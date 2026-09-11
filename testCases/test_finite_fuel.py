"""
# Finite-fuel software contracts

Compile the actual case for bounded initialisation, scalar-operator and
checkpoint probes. A separate four-step integration checks real restart
continuation. Reaction roots are compared to an independent bisection
evaluator of the implicit equation; none of these is a scientific campaign.
"""
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
PROBE = r'''
static int check_contract(void) {
  event("defaults");
  event("init");
  event("resume_state");
  if (exit_status || !fuel_ready) return exit_status ? exit_status : 9;
  if (!Period.x || !Period.y || !fuel.inverse ||
      rho1 != reynolds*density_ratio || rho2 != reynolds ||
      mu1 != viscosity_ratio || mu2 != 1. || stokes ||
      fuel_pe != pe || fuel_da != damkohler) return 10;
  double mass = fuel_inventory(), exterior = 0.;
  foreach(reduction(+:exterior)) exterior += (1. - f[])*dv();
  if (!resume_run && fabs(mass - initial_fuel*exterior) > 1e-10) return 11;
  fuel_surface_extension();
  if (!resume_run) {
    double error = 0.;
    foreach(reduction(max:error))
      if (fabs(d[]) < 2.*Delta)
        error = max(error, fabs(surface_fuel[] - initial_fuel));
    if (error > 1e-12) return 12;
  }
  const char * mode = getenv("FUEL_TEST");
  if (!mode) return 13;
  if (!strcmp(mode, "extension")) {
    for (int cycle = 0; cycle < 3; cycle++) {
      unrefine(level > min_level);
      refine(y > 0. && fabs(d[]) < 1. && level < max_level);
      if (!fuel_recover()) return 32;
      fuel_surface_extension();
      double error = 0.;
      foreach(reduction(max:error))
        if (fabs(d[]) < 2.*Delta)
          error = max(error, fabs(surface_fuel[] - initial_fuel));
      if (error > 1e-12) {
        fprintf(stderr, "AMR extension error=%g cycle=%d\n", error, cycle);
        return 33;
      }
    }
    return 0;
  }
  if (!strcmp(mode, "checkpoint")) {
    if (!resume_run) {
      foreach() if (f[] == 1.) fuel[] = DBL_EPSILON/4.;
      if (!fuel_repair_roundoff()) return 34;
    }
    if (!resume_run && !fuel_consume(0.2)) return 14;
    if (!resume_run) { t = 0.125; iter = 7; }
    if (!fuel_metadata(true)) return 15;
    dump(file = restart_file);
    printf("ACCOUNT %.17g %.17g %.17g t=%.17g i=%d roundoff=%.17g\n", fuel_inventory(), fuel_consumed,
           fuel_initial_inventory, t, iter, fuel_roundoff_inventory);
    return 0;
  }
  if (!strcmp(mode, "legacy")) {
    dump(file = restart_file, list = {f, d, u});
    return 0;
  }
  if (!strcmp(mode, "reaction")) {
    const double quantities[] = {0., 1e-20, 0.01, 1.};
    const double fractions[] = {1e-12, 0.01, 1.};
    const double exposures[] = {0., 1e-14, 0.3, 1e6};
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 4; k++) {
          double q = quantities[i]*fractions[j], h = fractions[j];
          double result = fuel_after_uptake(q, h, exposures[k], 0.1);
          if (!isfinite(result) || result < 0. || result > q) return 16;
          printf("ROOT %.17g %.17g %.17g %.17g\n", q, h, exposures[k], result);
        }
    return 0;
  }
  if (!strcmp(mode, "roundoff")) {
    fraction(f, 1. - sq(x) - sq(y));
    foreach() fuel[] = (1. - f[])*(0.5 + 0.49*sin(2.*pi*x/L0));
    boundary({f, fuel});
    double initial = fuel_inventory();
    fuel_initial_inventory = initial;
    fuel_roundoff_inventory = 0.;
    foreach_face(x) uf.x[] = 1.;
    foreach_face(y) uf.y[] = 0.73;
    boundary({uf});
    dt = 0.49*L0/(1 << max_level);
    for (int step = 0; step < 100; step++) {
      event("vof");
      if (exit_status || !fuel_recover()) return 35;
      if (fabs(fuel_inventory() - initial) > 1024.*DBL_EPSILON*initial) return 36;
    }
    if (!(fuel_roundoff_inventory > 0.) ||
        fuel_roundoff_inventory > 4096.*DBL_EPSILON*initial) return 37;
    fprintf(stderr, "ROUNDOFF steps=100 transferred=%.17g budget=%.17g balance=%.17g\n",
            fuel_roundoff_inventory, 4096.*DBL_EPSILON*initial,
            fuel_inventory() - initial);
    return 0;
  }
  if (!strcmp(mode, "reject_residue")) {
    for (int violation = 0; violation < 3; violation++) {
      foreach() {
        f[] = violation == 2 ? 1. + 6e-9 : (violation == 1 ? 0. : 1.);
        fuel[] = violation == 2 ? 0. : (violation == 1 ? -1e-8 : 1e-8);
      }
      double original = fuel_inventory(), ledger = fuel_roundoff_inventory;
      if (fuel_repair_roundoff() || fuel_inventory() != original ||
          fuel_roundoff_inventory != ledger) return 38;
    }
    return 0;
  }
  if (!strcmp(mode, "aligned")) {
    unrefine(level > min_level);
    for (int adaptive = 0; adaptive < 2; adaptive++) {
      if (adaptive) refine(x < 0. && y > 0. && level < max_level);
      foreach() {
        f[] = x < 0. ? 1. : 0.;
        fuel[] = 1. - f[];
      }
      boundary({f, fuel});
      fuel_consumed = 0.;
      double initial = fuel_inventory(), step = 1e-7;
      if (!fuel_consume(step)) return 30;
      double measured_area = fuel_consumed/(step*fuel_da/fuel_pe/(1. + fuel_half_saturation));
      if (fabs(measured_area - 2.*L0) > 1e-4 ||
          fabs(fuel_inventory() + fuel_consumed - initial) > 1e-10) return 31;
      printf("ALIGNED area=%.17g exact=%.17g adaptive=%d\n", measured_area, 2.*L0, adaptive);
    }
    return 0;
  }
  if (!strcmp(mode, "zero")) {
    if (initial_fuel == 0.) {
      refine(level < max_level);
      fraction(f, 1. - sq(x) - sq(y));
      foreach_face(x) uf.x[] = 1.;
      foreach_face(y) uf.y[] = 0.73;
      boundary({uf});
      dt = 0.49*L0/(1 << max_level);
      for (int step = 0; step < 5; step++) event("vof");
      if (exit_status || fuel_roundoff_inventory != 0.) return 39;
    }
    if (!fuel_diffuse(0.1) || !fuel_consume(10.)) return 17;
    if (initial_fuel == 0. && (fuel_inventory() != 0. || fuel_consumed != 0.)) return 18;
    if (damkohler == 0. && (fabs(fuel_inventory() - mass) > 1e-10 || fuel_consumed != 0.)) return 19;
    event("properties");
    double error = 0.;
    foreach(reduction(max:error))
      if (fabs(d[]) < 2.*Delta)
        error = max(error, fabs(sigmaf[] - (1./ca + gamma_slope*(1. - initial_fuel))));
    return error > 1e-10 ? 20 : 0;
  }
  if (!fuel_diffuse(0.1)) return 21;
  if (fabs(fuel_inventory() - mass) > 1e-10) return 22;
  foreach()
    fuel[] = (1. - f[])*(0.6 + 0.2*cos(2.*pi*x/L0) + 0.1*sin(2.*pi*y/L0));
  mass = fuel_inventory();
  for (int step = 0; step < 3; step++)
    if (!fuel_diffuse(0.03)) return 23;
  if (fabs(fuel_inventory() - mass) > 1e-8) return 24;
  foreach_face(x) uf.x[] = 0.2;
  foreach_face(y) uf.y[] = -0.1;
  boundary({uf});
  dt = 0.1*L0/(1 << max_level);
  for (int step = 0; step < 3; step++) vof_advection(interfaces, step);
  if (!fuel_recover() || fabs(fuel_inventory() - mass) > 1e-8) return 25;
  unrefine(level > min_level + 1);
  refine(f[] > 0. && f[] < 1. && level < max_level);
  if (!fuel_recover() || fabs(fuel_inventory() - mass) > 1e-8) return 26;
  event("adapt");
  if (exit_status || !fuel_recover() || fabs(fuel_inventory() - mass) > 1e-8) return 27;
  fuel_consumed = 0.;
  if (!fuel_consume(10.)) return 28;
  if (fuel_consumed <= 0. || fabs(fuel_inventory() + fuel_consumed - mass) > 1e-8) return 29;
  printf("PASS inventory=%.17g consumed=%.17g balance=%.17g\n",
         fuel_inventory(), fuel_consumed, fuel_inventory() + fuel_consumed - mass);
  return 0;
}
'''


class FiniteFuelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="finite-fuel-contract-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.work = Path(cls.temporary.name)
        compiler = REPO / "basilisk/src/qcc"
        if not compiler.is_file() or not (REPO / "basilisk/.comphy-lock").is_file():
            raise RuntimeError("Pinned project-local qcc is required")
        cls.environment = os.environ.copy()
        cls.environment["BASILISK"] = str(compiler.parent)
        cls.environment["PATH"] = str(compiler.parent) + os.pathsep + cls.environment.get("PATH", "")
        source = (REPO / "dropMove.c").read_text()
        source, count = re.subn(r"\brun\s*\(\s*\)\s*;", "return check_contract();", source)
        if count != 1:
            raise AssertionError("Expected one time-loop entry")
        source = source.replace("int main (", "static int check_contract(void);\nint main (", 1)
        (cls.work / "probe.c").write_text(source + PROBE)
        cls.executable = cls.work / "probe"
        compiled = subprocess.run([str(compiler), "-O1", "-disable-dimensions",
                                   f"-I{REPO}", "probe.c", "-o", str(cls.executable), "-lm"],
                                  cwd=cls.work, env=cls.environment, capture_output=True,
                                  text=True, timeout=120)
        if compiled.returncode:
            raise AssertionError(compiled.stdout + compiled.stderr)
        (cls.work / "actual.c").write_text((REPO / "dropMove.c").read_text())
        cls.actual = cls.work / "actual"
        compiled = subprocess.run([str(compiler), "-O1", "-disable-dimensions",
                                   f"-I{REPO}", "actual.c", "-o", str(cls.actual), "-lm"],
                                  cwd=cls.work, env=cls.environment, capture_output=True,
                                  text=True, timeout=120)
        if compiled.returncode:
            raise AssertionError(compiled.stdout + compiled.stderr)

    def invoke(self, mode="operators", parameters=(), success=True, directory=None):
        work = directory or Path(tempfile.mkdtemp(dir=self.work))
        environment = dict(self.environment, FUEL_TEST=mode)
        command = [str(self.executable), "min_level=3", "max_level=5", "domain=8",
                   "tmax=1", f"output_dir={work / 'output'}", *parameters]
        result = subprocess.run(command, cwd=work, env=environment,
                                capture_output=True, text=True, timeout=30)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_operators_inventory_and_phase_support(self):
        result = self.invoke()
        self.assertIn("PASS inventory=", result.stdout)

    def test_zero_fuel_and_zero_uptake(self):
        self.invoke("zero", ["InitialFuel=0"])
        self.invoke("zero", ["Da=0"])
        self.invoke("zero", ["Da=0", "GammaSlope=0", "InitialFuel=0.4"])

    def test_supplied_layouts_use_new_parameters(self):
        for name in ("single-drop", "two-drop", "four-drop"):
            with self.subTest(layout=name):
                self.invoke("zero", ["--params", str(REPO / "parameters" / f"{name}.cfg"),
                                     "domain=8", "max_level=5", "min_level=3", "Da=0",
                                     f"output_dir={self.work / name}"])

    def test_reaction_against_independent_bisection(self):
        result = self.invoke("reaction")
        rows = [row for row in result.stdout.splitlines() if row.startswith("ROOT ")]
        self.assertEqual(len(rows), 48)
        for row in rows:
            q, h, exposure, actual = map(float, row.split()[1:])
            lo, hi = 0., q
            for _ in range(180):
                trial = (lo + hi)/2.
                residual = trial + exposure*trial/(h*0.1 + trial) - q
                if residual > 0.:
                    hi = trial
                else:
                    lo = trial
            expected = (lo + hi)/2.
            self.assertAlmostEqual(actual, expected, delta=max(1e-300, q*2e-14))

    def test_face_aligned_interface_consumes_fuel_once(self):
        result = self.invoke("aligned")
        self.assertIn("ALIGNED area=", result.stdout)

    def test_constant_surface_extension_after_adaptation(self):
        self.invoke("extension", ["InitialFuel=0.7"])

    def test_high_cfl_vof_roundoff_is_conservative_and_accounted(self):
        result = self.invoke("roundoff", ["min_level=5", "max_level=5"])
        self.assertIn("ROUNDOFF steps=100", result.stderr)
        print(next(line for line in result.stderr.splitlines() if line.startswith("ROUNDOFF ")))
        self.invoke("reject_residue")

    def test_actual_timeloop_restart_advances_every_physical_step(self):
        def execute(directory, end, minimum, resume=False):
            command = [str(self.actual), "domain=8", f"min_level={minimum}", "max_level=5",
                       "tsnap=0.0001", f"tmax={end}", f"output_dir={directory}"]
            if resume:
                command.append("resume=1")
            result = subprocess.run(command, cwd=directory, env=self.environment,
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            rows = [list(map(float, line.split())) for line in
                    (directory / "log.dat").read_text().splitlines()[1:]]
            return rows[-1]
        for minimum in (5, 3):
            with self.subTest(min_level=minimum):
                full = Path(tempfile.mkdtemp(dir=self.work))
                split = Path(tempfile.mkdtemp(dir=self.work))
                reference = execute(full, "0.0004", minimum)
                checkpoint = execute(split, "0.0001", minimum)
                execute(split, "0.0002", minimum, True)
                resumed = execute(split, "0.0004", minimum, True)
                print("RESTART", "min_level=", minimum, "reference=", reference, "resumed=", resumed)
                self.assertAlmostEqual(resumed[1], reference[1], delta=1e-14)
                self.assertGreater(resumed[4], checkpoint[4]*3.9)
                self.assertAlmostEqual(resumed[4], reference[4], delta=1e-10)
                self.assertAlmostEqual(resumed[3], reference[3], delta=1e-9)
                self.assertAlmostEqual(resumed[2], reference[2],
                                       delta=max(1e-12, reference[2]*1e-6))

    def test_mapping_and_parameter_rejection(self):
        self.invoke("zero", ["Da=0", "Re=0.02", "densityRatio=2", "viscosityRatio=3"])
        for parameter in ("Oh=1", "AcNum=1", "Re=nan", "Da=-1", "Pe=0",
                          "FuelHalfSaturation=0", "InitialFuel=1.1", "Re=1junk",
                          "max_level=5.2", "GammaSlope=-1"):
            with self.subTest(parameter=parameter):
                self.invoke(parameters=[parameter], success=False)

    def test_checkpoint_preserves_inventory_and_rejects_legacy(self):
        work = Path(tempfile.mkdtemp(dir=self.work))
        original = self.invoke("checkpoint", directory=work)
        resumed = self.invoke("checkpoint", ["resume=1"], directory=work)
        self.assertEqual(original.stdout, resumed.stdout)
        self.invoke("checkpoint", ["resume=1", "Da=2"], directory=work, success=False)
        legacy = Path(tempfile.mkdtemp(dir=self.work))
        self.invoke("legacy", directory=legacy)
        self.invoke("checkpoint", ["resume=1"], directory=legacy, success=False)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""
Energy minimisation of a predicted structure with OpenMM.

Diffusion-based structure predictors produce near-ideal backbone torsions but
imperfect side-chain packing, which shows up in MolProbity as a high clashscore
alongside an implausibly good Ramachandran distribution. Restrained minimisation
relieves the steric overlaps while holding the backbone close to the predicted
conformation.

The backbone is restrained by default so that minimisation cannot quietly change
the fold that the rest of the analysis was performed on. Use --free to minimise
without restraints, but then re-check the superposition against the crystal
before reporting anything.

Usage:
    conda activate openmm
    python minimise.py input.pdb [-o output.pdb] [options]

Reports the RMSD from the starting coordinates so that the size of the change
is on record.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

try:
    from openmm import app, unit, LangevinMiddleIntegrator, CustomExternalForce, Platform
except ImportError:
    sys.exit("OpenMM not found. Run:  conda activate openmm")

try:
    from pdbfixer import PDBFixer
except ImportError:
    sys.exit("PDBFixer not found. Run:  conda install -c conda-forge pdbfixer")


FORCE_FIELDS = {
    "amber14":  ("amber14-all.xml", "amber14/tip3pfb.xml"),
    "charmm36": ("charmm36.xml", "charmm36/water.xml"),
}


def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("-o", "--output", default=None,
                   help="output PDB (default: <input>_min.pdb)")
    p.add_argument("--ff", choices=list(FORCE_FIELDS), default="amber14",
                   help="force field (default amber14)")
    p.add_argument("--restraint", type=float, default=10.0,
                   help="backbone positional restraint, kcal/mol/A^2 "
                        "(default 10.0; 0 disables)")
    p.add_argument("--free", action="store_true",
                   help="no restraints; equivalent to --restraint 0")
    p.add_argument("--tolerance", type=float, default=10.0,
                   help="convergence tolerance, kJ/mol/nm (default 10.0)")
    p.add_argument("--max-iterations", type=int, default=0,
                   help="0 = minimise until converged (default)")
    p.add_argument("--solvent", action="store_true",
                   help="minimise in explicit water rather than vacuum. "
                        "Slower; rarely changes the result for this purpose")
    p.add_argument("--keep-hydrogens", action="store_true",
                   help="keep input hydrogens instead of rebuilding them")
    args = p.parse_args()

    if args.free:
        args.restraint = 0.0

    src = Path(args.input).expanduser().resolve()
    if not src.is_file():
        sys.exit(f"Not found: {src}")
    out = (Path(args.output).expanduser().resolve() if args.output
           else src.with_name(src.stem + "_min.pdb"))

    print(f"Input       : {src.name}")
    print(f"Output      : {out.name}")
    print(f"Force field : {args.ff}")
    print(f"Restraint   : "
          f"{'none (free minimisation)' if args.restraint == 0 else f'{args.restraint} kcal/mol/A^2 on backbone'}")
    print(f"Environment : {'explicit solvent' if args.solvent else 'vacuum'}\n")

    # ---- prepare: add missing atoms and hydrogens
    print("Preparing structure...")
    fixer = PDBFixer(filename=str(src))
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # do not build unresolved loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    if fixer.missingAtoms:
        print(f"  adding {sum(len(v) for v in fixer.missingAtoms.values())} "
              f"missing heavy atoms")
    fixer.addMissingAtoms()
    if not args.keep_hydrogens:
        fixer.addMissingHydrogens(7.4)
        print("  hydrogens rebuilt at pH 7.4")

    topology, positions = fixer.topology, fixer.positions
    n_res = sum(1 for _ in topology.residues())
    print(f"  {n_res} residues, {topology.getNumAtoms()} atoms, "
          f"{topology.getNumChains()} chains")

    # ---- system
    ff_files = FORCE_FIELDS[args.ff]
    forcefield = app.ForceField(*ff_files)

    if args.solvent:
        print("\nAdding solvent...")
        modeller = app.Modeller(topology, positions)
        modeller.addSolvent(forcefield, padding=1.0 * unit.nanometer,
                            ionicStrength=0.15 * unit.molar)
        topology, positions = modeller.topology, modeller.positions
        print(f"  {topology.getNumAtoms()} atoms after solvation")
        system = forcefield.createSystem(
            topology, nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds)
    else:
        system = forcefield.createSystem(
            topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)

    # ---- positional restraints on backbone
    if args.restraint > 0:
        k = args.restraint * unit.kilocalories_per_mole / unit.angstrom ** 2
        force = CustomExternalForce(
            "0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
        force.addGlobalParameter("k", k)
        for name in ("x0", "y0", "z0"):
            force.addPerParticleParameter(name)
        n = 0
        for atom in topology.atoms():
            if (atom.name in ("CA", "C", "N", "O")
                    and atom.residue.name not in ("HOH", "NA", "CL")):
                force.addParticle(atom.index, positions[atom.index])
                n += 1
        system.addForce(force)
        print(f"\n  {n} backbone atoms restrained")

    # ---- minimise
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = Platform.getPlatformByName("CPU")
    sim = app.Simulation(topology, system, integrator, platform)
    sim.context.setPositions(positions)

    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"\nInitial energy : {e0.value_in_unit(unit.kilojoule_per_mole):14.1f} kJ/mol")

    t0 = time.time()
    sim.minimizeEnergy(
        tolerance=args.tolerance * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=args.max_iterations)

    elapsed = time.time() - t0

    state = sim.context.getState(getPositions=True, getEnergy=True)
    e1 = state.getPotentialEnergy()
    print(f"Final energy   : {e1.value_in_unit(unit.kilojoule_per_mole):14.1f} kJ/mol")
    print(f"Change         : {(e1 - e0).value_in_unit(unit.kilojoule_per_mole):14.1f} kJ/mol")
    print(f"Elapsed        : {elapsed:.1f} s")

    # ---- how far did it move?
    before = np.array([[v.x, v.y, v.z] for v in positions]) * 10.0   # nm -> A
    after = np.array([[v.x, v.y, v.z]
                      for v in state.getPositions()]) * 10.0

    heavy = [a.index for a in topology.atoms()
             if a.element is not None and a.element.symbol != "H"
             and a.residue.name not in ("HOH", "NA", "CL")]
    bb = [a.index for a in topology.atoms()
          if a.name in ("CA", "C", "N", "O")
          and a.residue.name not in ("HOH", "NA", "CL")]
    ca = [a.index for a in topology.atoms()
          if a.name == "CA" and a.residue.name not in ("HOH", "NA", "CL")]

    print(f"\nDisplacement from the input coordinates:")
    for label, idx in (("all heavy atoms", heavy), ("backbone", bb), ("CA only", ca)):
        if idx:
            d = np.linalg.norm(after[idx] - before[idx], axis=1)
            print(f"  {label:16s} RMSD {rmsd(after[idx], before[idx]):5.3f} A   "
                  f"max {d.max():5.3f} A")

    # ---- write, protein only
    if args.solvent:
        keep = [a.index for a in topology.atoms()
                if a.residue.name not in ("HOH", "NA", "CL")]
        sub = app.Modeller(topology, state.getPositions())
        sub.delete([r for r in topology.residues()
                    if r.name in ("HOH", "NA", "CL")])
        with open(out, "w") as fh:
            app.PDBFile.writeFile(sub.topology, sub.positions, fh, keepIds=True)
    else:
        with open(out, "w") as fh:
            app.PDBFile.writeFile(topology, state.getPositions(), fh, keepIds=True)

    print(f"\nWritten: {out}")
    print("""
Next steps
  1. Run MolProbity on the minimised structure and report the clashscore
     before and after. The Ramachandran statistics should change little.
  2. Re-run the superposition against the crystal to confirm the fold has
     not shifted. With backbone restraints the CA RMSD above should be
     well under 0.5 A.
  3. Report the validation measurements (burial, register, RMSD to crystal)
     from the UNMINIMISED model, and state that minimisation was applied
     only for downstream use.
""")


if __name__ == "__main__":
    main()

"""Small LM branch and resume checks without equilibrium solves."""

import json
from types import SimpleNamespace

import numpy as np
import pytest

import lsq_opt as lo


class BranchProblem:
    def __init__(self, *, jump=False, reject_first=False):
        self.jump = jump
        self.reject_first = reject_first
        self.calls = []

    def evaluate(self, x, mode="r", *, structure=None):
        x = np.asarray(x, float)
        self.calls.append((mode, structure))
        anchor = float(x[0]) if structure is None else structure.anchor
        fresh = SimpleNamespace(anchor=anchor) if structure is None else structure
        # Re-finding shifts this toy branch in the opposite direction to its
        # derivative, making premature convergence on a refound F observable.
        value = 1 + anchor - 2 * (x[0] - anchor) if self.jump else x[0] - 2
        slope = -2.0 if self.jump else 1.0
        ev = lo.Evaluation(x=x.copy(), boundary=x.copy(), ok=True,
                           r=np.array([value]), J=np.array([[slope]]) if mode == "J" else None)
        ev._chain = fresh
        if mode == "r" and self.reject_first:
            ev.ok = False
            self.reject_first = False
        return ev


def test_lm_freezes_trials_and_logs_refind_jump():
    p = BranchProblem(jump=True, reject_first=True)
    seen = []
    ev, state, message = lo.levenberg_marquardt(
        p, [0.0], maxiter=2, on_eval=lambda e, kind: seen.append((e, kind)))
    assert message == "maxiter"
    assert state["k"] == 2
    fresh = [e for e, kind in seen if kind == "J"]
    trials = [s for kind, s in p.calls if kind == "r"]
    assert trials[0] is trials[1] is fresh[0].structure
    assert trials[2] is fresh[1].structure
    assert all(s is None for kind, s in p.calls if kind == "J")
    jump = ev.info["structure_jump"]
    previous_trial = seen[-2][0]
    assert jump["material"]
    assert jump["F_delta"] == pytest.approx(ev.F - previous_trial.F)
    assert jump["residual_max"] == pytest.approx(abs(ev.r[0] - previous_trial.r[0]))


@pytest.mark.parametrize("tolerances", [dict(ftol=1e-10, xtol=0.0),
                                         dict(ftol=0.0, xtol=1e-8)], ids=["ftol", "xtol"])
def test_lm_material_refind_jump_prevents_convergence(tolerances):
    class RefindingJump(BranchProblem):
        def evaluate(self, x, mode="r", *, structure=None):
            ev = super().evaluate(x, mode, structure=structure)
            anchor = ev.structure.anchor
            ev.r = np.array([1 + 1e10 * anchor - (ev.x[0] - anchor)])
            if mode == "J":
                ev.J = np.array([[-1.0]])
            return ev

    p, seen = RefindingJump(), []
    ev, state, message = lo.levenberg_marquardt(
        p, [0.0], lam0=1e12, maxiter=2, **tolerances,
        on_eval=lambda e, kind: seen.append((e, kind)))
    initial, trial, refound = [e for e, _ in seen[:3]]
    assert refound.info["structure_jump"]["material"]
    if tolerances["ftol"]:
        assert initial.F - trial.F <= tolerances["ftol"] * trial.F
    else:
        assert np.linalg.norm(trial.x - initial.x) <= tolerances["xtol"] * (1 + np.linalg.norm(trial.x))
    assert message == "maxiter"
    assert state["k"] == 2
    assert sum(kind == "J" for kind, _ in p.calls) == 3
    assert ev.info["structure_jump"]["material"]


def test_lm_resume_recreates_accepted_structure():
    kwargs = dict(maxiter=3, lam0=1.0, ftol=0.0, xtol=0.0)
    full, full_state, _ = lo.levenberg_marquardt(BranchProblem(), [0.0], **kwargs)
    _, saved, _ = lo.levenberg_marquardt(BranchProblem(), [0.0], **dict(kwargs, maxiter=1))
    saved = json.loads(json.dumps(saved, default=lambda a: np.asarray(a).tolist()))
    p = BranchProblem()
    resumed, resumed_state, _ = lo.levenberg_marquardt(p, [0.0], state=saved, **kwargs)
    assert p.calls[0] == ("J", None)
    np.testing.assert_allclose(resumed.x, full.x, rtol=0, atol=0)
    assert resumed_state["lam"] == full_state["lam"]
    assert resumed_state["k"] == full_state["k"]


def test_problem_refinds_after_frozen_trial_without_resolving(monkeypatch):
    class Chain:
        def __init__(self, *args, harm):
            self.anchor = float(harm["x"][0])

        def values(self, harm):
            return np.array([[1 + self.anchor - 2 * (float(harm["x"][0]) - self.anchor)]])

    monkeypatch.setattr(lo, "WoutChain", Chain)
    monkeypatch.setattr(lo, "term_values", lambda values, objective: values)
    monkeypatch.setattr(lo.vj.STATS, "last", {"converged": True})
    p = object.__new__(lo.Problem)
    p.base_wout, p.s, p.settings, p.mboz, p.nboz = None, None, None, None, None
    p.n_eval, p.objective = 0, "test"
    p.reference = SimpleNamespace(term_sums=np.ones(1))
    p.boundary = lambda x: np.asarray(x)
    p._harm = lambda boundary: {"x": boundary}
    p.params = lambda: {"mode": "pen"}
    p._residual_jit = lambda mode: (lambda values, *args: values.ravel(), None)
    p.constraints = lambda *args: (np.zeros(0), np.zeros(0))

    base = p.evaluate([0.0])
    trial = p.evaluate([0.5], structure=base.structure)
    assert trial.r[0] == pytest.approx(0.0)
    refound = p.evaluate([0.5])
    assert refound.r[0] == pytest.approx(1.5)
    assert refound.structure is not trial.structure
    assert p.n_eval == 2
    assert refound.info["reused_forward"]
    assert not refound.info["frozen_structure"]
    # Switching back to an older branch at the same boundary also invalidates
    # values, while retaining the same cached equilibrium.
    again = p.evaluate([0.5], structure=base.structure)
    assert again.r[0] == pytest.approx(0.0)
    assert p.n_eval == 2

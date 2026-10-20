"""Replay tier-changing PR events through the workflows' trigger filters."""

import fnmatch
from pathlib import Path
import unittest

import yaml


WORKFLOWS = Path(__file__).resolve().parents[1] / ".github" / "workflows"


def dispatches(workflow, action, base):
    # GitHub's pull_request dispatch defaults and branch-filter semantics:
    # https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#pull_request
    events = workflow["on"]
    if "pull_request" not in events:
        return False
    filters = events["pull_request"] or {}
    types = filters.get("types", ["opened", "synchronize", "reopened"])
    branches = filters.get("branches", ["**"])
    return action in types and any(fnmatch.fnmatchcase(base, b) for b in branches)


class FullTierEvents(unittest.TestCase):
    def check_dispatch(self, action, bases):
        for name in ("benchmark", "golden"):
            with (WORKFLOWS / f"{name}.yml").open() as stream:
                workflow = yaml.load(stream, Loader=yaml.BaseLoader)
            for base in bases:
                with self.subTest(workflow=name, action=action, base=base):
                    self.assertTrue(dispatches(workflow, action, base))

    def test_retarget_reconsiders_full_tier(self):
        # Retargeting to main must run; retargeting to a stack must report skip.
        self.check_dispatch("edited", ("main", "stack-base"))

    def test_label_removal_reconsiders_full_tier(self):
        # Removing full-ci must re-evaluate the job, even without a new commit.
        self.check_dispatch("unlabeled", ("main", "stack-base"))


if __name__ == "__main__":
    unittest.main()

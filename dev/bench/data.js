window.BENCHMARK_DATA = {
  "lastUpdate": 1779113687596,
  "repoUrl": "https://github.com/itpplasma/rabe",
  "entries": {
    "Benchmark": [
      {
        "commit": {
          "author": {
            "email": "100583655+GeorgGrassler@users.noreply.github.com",
            "name": "GeorgGrassler",
            "username": "GeorgGrassler"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cdbff9ad2a86ed770d6a855255bfd64cc8b7fa2c",
          "message": "remove: python (#65)",
          "timestamp": "2026-05-18T13:49:21+02:00",
          "tree_id": "d3ffb1f0492357e94c25a4d7bb49f855c6eb6d36",
          "url": "https://github.com/itpplasma/rabe/commit/cdbff9ad2a86ed770d6a855255bfd64cc8b7fa2c"
        },
        "date": 1779105217501,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 1.0021,
            "unit": "ratio",
            "extra": "current: 58160 ms, baseline: 58040 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 58160,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "100583655+GeorgGrassler@users.noreply.github.com",
            "name": "GeorgGrassler",
            "username": "GeorgGrassler"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "eedaf245c2704b92c92862476d009a2bbb73dd0b",
          "message": "Update/quadpack (#64)\n\n* use: modern quadpack\n\n* update: golden record\n\n* update: README\n\n* update: third party dir name",
          "timestamp": "2026-05-18T15:25:15+02:00",
          "tree_id": "a8ba4577dddf1d4a8636b1641b8b712dc193bc1f",
          "url": "https://github.com/itpplasma/rabe/commit/eedaf245c2704b92c92862476d009a2bbb73dd0b"
        },
        "date": 1779110972501,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 1.0126,
            "unit": "ratio",
            "extra": "current: 59570 ms, baseline: 58830 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 59570,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "100583655+GeorgGrassler@users.noreply.github.com",
            "name": "GeorgGrassler",
            "username": "GeorgGrassler"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "45717986d1bb5079e703930ddf9810ed1aeb7198",
          "message": "Fix/git in build (#66)\n\n* use: main in case of empty branch\n\nIf the current directory is not a git repo, BRANCH is empty but still\nleads to a BRANCH_EXISTS as git ls-remote --heads ${REPO_URL} ${BRANCH}\ngives back all remote branches without failure.\n\n* let: git fail gracefully\n\nIf the directory is not a git repo, the build now still goes through,\nbut sets \"unknown\" as the git hash.\n\n* add: hard git dependance\n\nThe build process requires git for fetching dependencies.",
          "timestamp": "2026-05-18T16:10:40+02:00",
          "tree_id": "14147fa9ab5417fcb4b8182331a8ed798d762483",
          "url": "https://github.com/itpplasma/rabe/commit/45717986d1bb5079e703930ddf9810ed1aeb7198"
        },
        "date": 1779113687026,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 1.004,
            "unit": "ratio",
            "extra": "current: 59530 ms, baseline: 59290 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 59530,
            "unit": "ms"
          }
        ]
      }
    ]
  }
}
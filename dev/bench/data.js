window.BENCHMARK_DATA = {
  "lastUpdate": 1781113187654,
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
          "id": "550ccebd7993311d30e2598db981bc0624dbff2e",
          "message": "update/doc (#69)\n\n* added plotting scripts to golden record example\n\n* extended README",
          "timestamp": "2026-05-19T14:52:44+02:00",
          "tree_id": "22e5263cde53b19828ef70a4dcbcb59f4e591f2e",
          "url": "https://github.com/itpplasma/rabe/commit/550ccebd7993311d30e2598db981bc0624dbff2e"
        },
        "date": 1779195397696,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 0.9962,
            "unit": "ratio",
            "extra": "current: 68150 ms, baseline: 68410 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 68150,
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
          "id": "92dce3845f77fd887df86246c92894c1609646c0",
          "message": "check: origin with tol (#67)\n\n* check: origin with tol\n\nDue to the field slightly violating stellator symmetry, the minimum\nmight not lie at the minimum (including error tolerances). We therefore\ncheck the the B value difference between the origin and the found\nminimum is below the symmetry violation.\n\n* archive: git snapeshot",
          "timestamp": "2026-05-19T17:19:10+02:00",
          "tree_id": "677f0c24b01393e3d23f42fe95ca8a59b23970c6",
          "url": "https://github.com/itpplasma/rabe/commit/92dce3845f77fd887df86246c92894c1609646c0"
        },
        "date": 1779204184312,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 0.9922,
            "unit": "ratio",
            "extra": "current: 71040 ms, baseline: 71600 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 71040,
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
          "id": "e9e2a47ebe446ab00c97c14d681f8a73de48eb82",
          "message": "Release/improvements (#70)\n\n* write: dat file and update README\n\n* specify: surfaces with range\n\n* rename: output coefficients\n\n* use: explicit env\n\nwe only use it for libneo and want to avoid overlap with other common env names\n\n* change: executable extension\n\n* update: README\n\n* make: build generator choice explicit\n\n* export: PDF on tag release\n\n* update: golden record and adjust tolerances",
          "timestamp": "2026-05-20T18:37:28+02:00",
          "tree_id": "eca7aae813ae9db07ed199ea16001042a4796a18",
          "url": "https://github.com/itpplasma/rabe/commit/e9e2a47ebe446ab00c97c14d681f8a73de48eb82"
        },
        "date": 1779295315160,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 1.0182,
            "unit": "ratio",
            "extra": "current: 72580 ms, baseline: 71280 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 72580,
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
          "id": "efcf0b27e60479fb1f1155648661122ae067eb43",
          "message": "specify: tolerance as field property (#71)\n\nDepending on the field type, the numerical value of B can deviate from the physics value of B. This puts a limit on how good certain physics properties like stellarator symmetry can be fulfilled. The information is bound to the field.",
          "timestamp": "2026-05-21T10:20:52+02:00",
          "tree_id": "ae2af2a7800c5b91b34a9c8199ed82276fc22d1a",
          "url": "https://github.com/itpplasma/rabe/commit/efcf0b27e60479fb1f1155648661122ae067eb43"
        },
        "date": 1779351898836,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 0.9943,
            "unit": "ratio",
            "extra": "current: 71530 ms, baseline: 71940 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 71530,
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
          "id": "a8c3b21d89dc46fb63faeb172d58f7b62a3dd206",
          "message": "58 prepare codebase for public release remove mandatory suitesparse dependency (#74)\n\n* replace: neo_field_t by analytic fields\n\nAnd correct analytical estimate formula.\n\n* move: neo_field_t to rabe_internal\n\n* move: data/field files to rabe_internal\n\n* use: explicit CMake dir locations\n\n- PROJECT_* instead of CMAKE_* i.e. always in respect to the project CMakeList\n- CMAKE_CURRENT_LIST_DIR for relative to the current CMakeList\n- CMAKE_CURRENT_BINARY_DIR for the mirrored project location in the build directory (i.e. the patched .git_version)\n\n* add: internal as submodul\n\n* update: docs",
          "timestamp": "2026-05-22T17:32:17+02:00",
          "tree_id": "327eab7496399af75c26be24ac6b6eb86375661f",
          "url": "https://github.com/itpplasma/rabe/commit/a8c3b21d89dc46fb63faeb172d58f7b62a3dd206"
        },
        "date": 1779464186478,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 0.9905,
            "unit": "ratio",
            "extra": "current: 71170 ms, baseline: 71850 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 71170,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "albert@tugraz.at",
            "name": "Christopher Albert",
            "username": "krystophny"
          },
          "committer": {
            "email": "albert@tugraz.at",
            "name": "Christopher Albert",
            "username": "krystophny"
          },
          "distinct": true,
          "id": "2fac90c7d2e55aad911afbceab0249649fe4774c",
          "message": "Widen golden-record tolerance for cross-platform runs\n\nmacOS ARM and Linux x86 builds diverge at ~3e-9 relative for most\noutputs and ~4e-7 for Lambda_A through libm and adaptive-quadrature\ndifferences. rtol=1e-10 failed on any machine other than the one\nthat produced the expected record; 1e-6 still catches algorithmic\nregressions, which start at 1e-3 in practice.",
          "timestamp": "2026-06-10T19:36:45+02:00",
          "tree_id": "053728382e77ea0b9c7824f95377001d6932b117",
          "url": "https://github.com/itpplasma/rabe/commit/2fac90c7d2e55aad911afbceab0249649fe4774c"
        },
        "date": 1781113186698,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 1.0005,
            "unit": "ratio",
            "extra": "current: 38580 ms, baseline: 38560 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 38580,
            "unit": "ms"
          }
        ]
      }
    ]
  }
}
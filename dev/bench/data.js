window.BENCHMARK_DATA = {
  "lastUpdate": 1795897879345,
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
      },
      {
        "commit": {
          "author": {
            "email": "albert@tugraz.at",
            "name": "Christopher Albert",
            "username": "krystophny"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "baabd080e77a0a46d8630ae929fff3b297dba336",
          "message": "Revert quadpack link and golden tolerance changes (#86)\n\nReverts 32566ea and 2fac90c.",
          "timestamp": "2026-06-10T20:09:36+02:00",
          "tree_id": "18233626bf914324cb5dd1adc4de5a51fccd810e",
          "url": "https://github.com/itpplasma/rabe/commit/baabd080e77a0a46d8630ae929fff3b297dba336"
        },
        "date": 1781115213754,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 0.9956,
            "unit": "ratio",
            "extra": "current: 70650 ms, baseline: 70960 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 70650,
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
          "id": "8d14bbc6506b28af02ac455bed755abde93bb18f",
          "message": "Merge pull request #85 from GeorgGrassler/main\n\nupdate: upstream",
          "timestamp": "2026-06-10T20:32:35+02:00",
          "tree_id": "821fc7253b64c76b7dccabce8d2c7894c22a3ec6",
          "url": "https://github.com/itpplasma/rabe/commit/8d14bbc6506b28af02ac455bed755abde93bb18f"
        },
        "date": 1781116605200,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 1.0027,
            "unit": "ratio",
            "extra": "current: 71470 ms, baseline: 71280 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 71470,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "georg.grassler+bot@student.tugraz.at",
            "name": "gg-codebot",
            "username": "gg-codebot"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": false,
          "id": "94d849409e5927afc13f94c483aa294f30901b99",
          "message": "ci: run checks in the merge queue, Benchmark history on main again (#149)\n\n* ci: run checks in the merge queue, Benchmark history on main again\n\ntest, Golden Record and Benchmark run on merge_group. Benchmark runs in\nthe queue and on push to main (stores the history again, lost when the\npush trigger was removed); on PRs its job is skipped.\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n\n* ci: Benchmark runs on PRs again\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n\n* Require the maintainer's review via CODEOWNERS\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n\n---------\n\nCo-authored-by: Bot <bot@Testers-MacBook-Pro.local>\nCo-authored-by: Claude Opus 5.5 <noreply@anthropic.com>",
          "timestamp": "2026-11-27T18:32:46Z",
          "tree_id": "b5e23345ac2b336b1049d67e7f61022eced90d66",
          "url": "https://github.com/itpplasma/rabe/commit/94d849409e5927afc13f94c483aa294f30901b99"
        },
        "date": 1795805114175,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 0.9968,
            "unit": "ratio",
            "extra": "current: 31060 ms, baseline: 31160 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 31060,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "georg.grassler+bot@student.tugraz.at",
            "name": "gg-codebot",
            "username": "gg-codebot"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": false,
          "id": "2f375d0c5d812005e84b1242ec7cac4798146363",
          "message": "feat: boozmn spline grid as input; electric_rabe loader by field_type (#155)\n\n* feat: boozmn spline grid as input (n_rho, n_theta, n_zeta)\n\nThe boozmn loader used libneo's fixed 30 x 48 x 96 grid. For CIEMAT QI,\nderivatives of |B| change by 0.6 % when n_theta goes from 48 to 96.\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n(cherry picked from commit 5ea10066a0b97c8f22302a7a3d8353e41332ea3f)\n\n* feat: electric_rabe selects the field loader by field_type\n\nSame selection as rabe.x, so electric_rabe can read boozmn files.\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n(cherry picked from commit a10d4a0ee73e2e7b3197aabbcd3383ea3b0bd508)\n\n* update: Readme\n\n---------\n\nCo-authored-by: Bot <bot@nt-b070.wlcli.tugraz.at>\nCo-authored-by: Claude Opus 5.5 <noreply@anthropic.com>\nCo-authored-by: GeorgGrassler <100583655+GeorgGrassler@users.noreply.github.com>",
          "timestamp": "2026-11-28T18:00:12Z",
          "tree_id": "0fcfcbe767962de1f19cf30203eee646da1aa2f7",
          "url": "https://github.com/itpplasma/rabe/commit/2f375d0c5d812005e84b1242ec7cac4798146363"
        },
        "date": 1795889524067,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 0.998,
            "unit": "ratio",
            "extra": "current: 25460 ms, baseline: 25510 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 25460,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "georg.grassler+bot@student.tugraz.at",
            "name": "gg-codebot",
            "username": "gg-codebot"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": false,
          "id": "0f8812af9200b0b14acc545ad4f2f0b5aa5ef55c",
          "message": "fix: electric_rabe applies the covariant factor twice in prefactor_A (on main) (#151)\n\n* fix: electric_rabe applies the covariant factor twice in prefactor_A\n\nflock%I_ref contains the covariant factor since 857dc93 (Handle covariant\ncomponents); get_g_modes_from_fieldlines multiplied by it again, so the\nboundary-layer term (1/sqrt(nu_star)) was too large by sqrt(B_phi + iota B_theta).\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n\n* remove: comment\n\n---------\n\nCo-authored-by: bot <bot@tugraz.at>\nCo-authored-by: Claude Opus 5.5 <noreply@anthropic.com>\nCo-authored-by: Bot <bot@Testers-MacBook-Pro.local>\nCo-authored-by: GeorgGrassler <100583655+GeorgGrassler@users.noreply.github.com>",
          "timestamp": "2026-11-28T18:25:05Z",
          "tree_id": "8ec1a532a0939b3f0197846b78f099745de9bb64",
          "url": "https://github.com/itpplasma/rabe/commit/0f8812af9200b0b14acc545ad4f2f0b5aa5ef55c"
        },
        "date": 1795890955817,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 0.9916,
            "unit": "ratio",
            "extra": "current: 30840 ms, baseline: 31100 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 30840,
            "unit": "ms"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "georg.grassler+bot@student.tugraz.at",
            "name": "gg-codebot",
            "username": "gg-codebot"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": false,
          "id": "21434fc92d1a29429aab01fe2102f9c7258c7be3",
          "message": "Optional field line diagnostic per surface: flocks_dir/flock_<k>.nc (#152)\n\n* feat: rabe writes the field lines per surface to fieldlines.nc\n\nLabels, maxima bounding the well and delta_eta of each field line, so the\nviolation of omnigenity at the maxima can be inspected.\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n\n* feat: field line diagnostic per surface into flocks_dir instead of fieldlines.nc\n\nIf flocks_dir is set, rabe checks at startup that it is writable and writes\neach flock to flocks_dir/flock_<k>.nc right after it is computed: the flock\nscalars and the field line components along their own dimension, no padding.\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n\n* fix: missing continuation after flocks_dir in rabe_config from the main merge\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n\n* udpate: docs\n\n* add: dir as input\n\n---------\n\nCo-authored-by: Bot <bot@nt-b240.wlcli.tugraz.at>\nCo-authored-by: Claude Opus 5.5 <noreply@anthropic.com>\nCo-authored-by: Bot <bot@nt-f144.wlcli.tugraz.at>\nCo-authored-by: Bot <bot@Testers-MacBook-Pro.local>\nCo-authored-by: GeorgGrassler <100583655+GeorgGrassler@users.noreply.github.com>\nCo-authored-by: Georg Grassler <georg.grassler@tugraz.at>",
          "timestamp": "2026-11-28T20:20:27Z",
          "tree_id": "46ca5e349ee80a26aebe2578da79a95060c9fe99",
          "url": "https://github.com/itpplasma/rabe/commit/21434fc92d1a29429aab01fe2102f9c7258c7be3"
        },
        "date": 1795897878548,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "execution time ratio (current/baseline)",
            "value": 0.9892,
            "unit": "ratio",
            "extra": "current: 15570 ms, baseline: 15740 ms"
          },
          {
            "name": "execution time (current, Release)",
            "value": 15570,
            "unit": "ms"
          }
        ]
      }
    ]
  }
}
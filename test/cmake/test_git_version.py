"""Compile and execute embedded metadata after ordinary Git/ref updates."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile


parser = argparse.ArgumentParser()
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--cmake', required=True)
parser.add_argument('--git', required=True)
parser.add_argument('--compiler', required=True)
parser.add_argument('--generator', required=True)
args = parser.parse_args()


def run(command, cwd=None):
    return subprocess.check_output(command, cwd=cwd, text=True,
                                   stderr=subprocess.STDOUT).strip()


def git(root, *arguments):
    return run([args.git, '-C', str(root), *arguments])


def build_and_check(source, build):
    run([args.cmake, '--build', str(build), '--target', 'version_probe'])
    expected = git(source, 'rev-parse', 'HEAD')
    actual = run([str(build / 'version_probe')])
    assert actual == expected, (source, expected, actual)
    return actual


def configure(source, build):
    run([args.cmake, '-S', str(source), '-B', str(build), '-G', args.generator,
         '-DCMAKE_Fortran_COMPILER=' + args.compiler])
    return build_and_check(source, build)


def commit(source, message):
    # This repository contains only a synthetic CMake/Fortran fixture. Use the
    # already configured identity where present, with a fixture identity for CI.
    try:
        name = run([args.git, 'config', '--get', 'user.name'], args.source)
        email = run([args.git, 'config', '--get', 'user.email'], args.source)
    except subprocess.CalledProcessError:
        name, email = 'RABE regression fixture', 'fixture@example.invalid'
    git(source, 'config', 'user.name', name)
    git(source, 'config', 'user.email', email)
    git(source, 'add', 'CMakeLists.txt', 'main.f90', 'src/git_version.f90.in',
        'cmake/git_version.cmake', 'revision')
    git(source, 'commit', '-m', message)
    return git(source, 'rev-parse', 'HEAD')


with tempfile.TemporaryDirectory(prefix='rabe-git-version-') as temporary:
    base = Path(temporary)
    source = base / 'main'
    source.mkdir()
    git(source, 'init', '-b', 'main')
    (source / 'src').mkdir()
    (source / 'cmake').mkdir()
    shutil.copyfile(args.source / 'cmake/git_version.cmake',
                    source / 'cmake/git_version.cmake')
    shutil.copyfile(args.source / 'src/git_version.f90.in',
                    source / 'src/git_version.f90.in')
    module = (args.source / 'cmake/git_version_target.cmake').as_posix()
    (source / 'CMakeLists.txt').write_text(f'''cmake_minimum_required(VERSION 3.24)
project(version_fixture LANGUAGES Fortran)
include("{module}")
add_library(version_lib "${{PROJECT_BINARY_DIR}}/src/git_version.f90")
add_dependencies(version_lib git_version_target)
add_executable(version_probe main.f90)
target_link_libraries(version_probe PRIVATE version_lib)
''')
    (source / 'main.f90').write_text('''program version_probe
    use git_version, only: git_hash
    implicit none
    print '(a)', git_hash
end program version_probe
''')
    (source / 'revision').write_text('initial\n')
    first = commit(source, 'initial fixture')
    main_build = base / 'main-build'
    configure(source, main_build)
    (source / 'revision').write_text('second\n')
    second = commit(source, 'advance regular checkout')
    assert first != second
    build_and_check(source, main_build)

    # Moving the active branch ref does not rewrite .git/HEAD or the index.
    git(source, 'update-ref', 'refs/heads/main', first)
    build_and_check(source, main_build)
    git(source, 'update-ref', 'refs/heads/main', second)
    build_and_check(source, main_build)
    print('regular checkout: commit and active-ref changes refresh embedded hash')

    linked = base / 'linked'
    git(source, 'worktree', 'add', '-b', 'linked', str(linked))
    try:
        linked_build = base / 'linked-build'
        configure(linked, linked_build)
        (linked / 'revision').write_text('linked\n')
        third = commit(linked, 'advance linked worktree')
        build_and_check(linked, linked_build)
        git(source, 'pack-refs', '--all')
        git(linked, 'update-ref', 'refs/heads/linked', first)
        build_and_check(linked, linked_build)
        git(linked, 'update-ref', 'refs/heads/linked', third)
        build_and_check(linked, linked_build)
        git(linked, 'checkout', '--detach', second)
        build_and_check(linked, linked_build)
        print('linked worktree: commit, packed ref and detached HEAD refresh hash')

        # A no-op revision check must leave generated source and executable
        # untouched, so incremental builds do not recompile or relink them.
        generated = linked_build / 'src/git_version.f90'
        executable = linked_build / 'version_probe'
        before = (generated.stat().st_mtime_ns, executable.stat().st_mtime_ns)
        build_and_check(linked, linked_build)
        after = (generated.stat().st_mtime_ns, executable.stat().st_mtime_ns)
        assert before == after, ('no-op build rewrote metadata or relinked', before, after)
        print('unchanged HEAD: generated source and executable timestamps unchanged')
    finally:
        git(source, 'worktree', 'remove', '--force', str(linked))

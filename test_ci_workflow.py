#!/usr/bin/env python3
"""
Test script to verify GitHub Actions CI workflow configuration.
This is a RED phase test - it should fail until the workflow is implemented.
"""

import os
import sys
import yaml
from pathlib import Path


def test_workflow_file_exists():
    """Test that the CI workflow file exists in the correct location."""
    workflow_path = Path(".github/workflows/ci.yml")
    assert workflow_path.exists(), f"Workflow file not found at {workflow_path}"
    return workflow_path


def test_workflow_structure(workflow_path):
    """Test that the workflow has the correct structure."""
    with open(workflow_path, 'r') as f:
        workflow = yaml.safe_load(f)
    
    # Check workflow name
    assert 'name' in workflow, "Workflow must have a name"
    assert workflow['name'] == 'CI', f"Workflow name should be 'CI', got '{workflow['name']}'"
    
    # Check triggers
    assert 'on' in workflow, "Workflow must have triggers"
    triggers = workflow['on']
    
    # Check push trigger
    assert 'push' in triggers, "Workflow must trigger on push"
    # Push should trigger on all branches (no branch filter)
    
    # Check pull_request trigger
    assert 'pull_request' in triggers, "Workflow must trigger on pull_request"
    pr_config = triggers['pull_request']
    assert 'branches' in pr_config, "Pull request must specify target branches"
    assert 'main' in pr_config['branches'], "Pull request must target main branch"
    
    # Check jobs
    assert 'jobs' in workflow, "Workflow must have jobs"
    return workflow


def test_build_matrix(workflow):
    """Test that the workflow has the correct build matrix."""
    jobs = workflow['jobs']
    
    # Check for test job
    assert 'test' in jobs, "Workflow must have a 'test' job"
    test_job = jobs['test']
    
    # Check strategy
    assert 'strategy' in test_job, "Test job must have a strategy"
    strategy = test_job['strategy']
    
    # Check matrix
    assert 'matrix' in strategy, "Strategy must have a matrix"
    matrix = strategy['matrix']
    
    # Check OS matrix
    assert 'os' in matrix, "Matrix must include OS"
    os_list = matrix['os']
    assert 'ubuntu-latest' in os_list, "Matrix must include ubuntu-latest"
    assert 'macos-latest' in os_list, "Matrix must include macos-latest"
    
    # Check build type matrix
    assert 'build-type' in matrix, "Matrix must include build-type"
    build_types = matrix['build-type']
    assert 'Debug' in build_types, "Matrix must include Debug build"
    assert 'Release' in build_types, "Matrix must include Release build"
    
    return test_job


def test_job_steps(test_job):
    """Test that the job has all required steps."""
    assert 'steps' in test_job, "Test job must have steps"
    steps = test_job['steps']
    
    # Check for essential steps
    step_names = [step.get('name', '') for step in steps]
    
    required_steps = [
        'Checkout code',
        'Set up Python',
        'Install system dependencies',
        'Install Python dependencies',
        'Configure CMake',
        'Build',
        'Run quick tests',
    ]
    
    for required in required_steps:
        assert any(required in name for name in step_names), f"Missing required step: {required}"


def test_slow_test_job(workflow):
    """Test that there's a separate job for slow tests."""
    jobs = workflow['jobs']
    
    # Check for slow test job
    assert 'test-slow' in jobs, "Workflow must have a 'test-slow' job"
    slow_job = jobs['test-slow']
    
    # Check it only runs on pull requests
    assert 'if' in slow_job, "Slow test job must have a condition"
    condition = slow_job['if']
    assert 'pull_request' in condition, "Slow tests should only run on pull requests"
    
    # Check it only runs on ubuntu-latest
    if 'strategy' in slow_job and 'matrix' in slow_job['strategy']:
        matrix = slow_job['strategy']['matrix']
        if 'os' in matrix:
            assert matrix['os'] == ['ubuntu-latest'], "Slow tests should only run on ubuntu-latest"
    
    # Check for timeout
    assert 'timeout-minutes' in slow_job or any(
        'timeout-minutes' in step for step in slow_job.get('steps', [])
    ), "Slow test job should have a timeout"


def test_python_tests_integration(test_job):
    """Test that Python tests are integrated."""
    steps = test_job['steps']
    
    # Check for Python test execution
    step_commands = []
    for step in steps:
        if 'run' in step:
            step_commands.append(step['run'])
    
    assert any('pytest' in cmd or 'python -m pytest' in cmd for cmd in step_commands), \
        "Workflow must run Python tests with pytest"


def test_ctest_label_filtering(test_job):
    """Test that CTest is configured with proper label filtering."""
    steps = test_job['steps']
    
    # Check for CTest execution with labels
    for step in steps:
        if 'name' in step and 'quick test' in step['name'].lower():
            assert 'run' in step, "Quick test step must have a run command"
            run_cmd = step['run']
            assert 'ctest' in run_cmd, "Must use ctest for running tests"
            assert '-L quick' in run_cmd or '--label-regex quick' in run_cmd, \
                "Quick tests must be filtered by label"


def main():
    """Run all tests."""
    print("Running CI workflow tests...")
    failures = []
    
    try:
        # Test 1: File exists
        print("✓ Testing workflow file existence...")
        workflow_path = test_workflow_file_exists()
        print("  PASS: Workflow file exists")
        
        # Test 2: Workflow structure
        print("✓ Testing workflow structure...")
        workflow = test_workflow_structure(workflow_path)
        print("  PASS: Workflow structure is correct")
        
        # Test 3: Build matrix
        print("✓ Testing build matrix...")
        test_job = test_build_matrix(workflow)
        print("  PASS: Build matrix is configured correctly")
        
        # Test 4: Job steps
        print("✓ Testing job steps...")
        test_job_steps(test_job)
        print("  PASS: All required steps are present")
        
        # Test 5: Slow test job
        print("✓ Testing slow test job...")
        test_slow_test_job(workflow)
        print("  PASS: Slow test job is configured correctly")
        
        # Test 6: Python tests
        print("✓ Testing Python test integration...")
        test_python_tests_integration(test_job)
        print("  PASS: Python tests are integrated")
        
        # Test 7: CTest labels
        print("✓ Testing CTest label filtering...")
        test_ctest_label_filtering(test_job)
        print("  PASS: CTest label filtering is configured")
        
    except AssertionError as e:
        print(f"  FAIL: {e}")
        failures.append(str(e))
    except Exception as e:
        print(f"  ERROR: {e}")
        failures.append(str(e))
    
    if failures:
        print(f"\n❌ {len(failures)} test(s) failed:")
        for failure in failures:
            print(f"  - {failure}")
        sys.exit(1)
    else:
        print("\n✅ All tests passed!")
        sys.exit(0)


if __name__ == "__main__":
    main()
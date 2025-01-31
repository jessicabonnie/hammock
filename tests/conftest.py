import pytest # type: ignore

def pytest_configure(config):
    """Add markers to the pytest configuration."""
    config.addinivalue_line("markers", "quick: mark test as quick to run")
    config.addinivalue_line("markers", "full: mark test as part of the full test suite")
    config.addinivalue_line("markers", "slow: mark test as very slow to run") 
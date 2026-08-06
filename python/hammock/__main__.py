import sys

from hammock.cli import main

if __name__ == "__main__":
    # sys.exit(...) matters: main() returns the exit status, and without it
    # `python -m hammock` reports success on every error path (the console
    # script installed from pyproject already wraps main() this way).
    sys.exit(main())

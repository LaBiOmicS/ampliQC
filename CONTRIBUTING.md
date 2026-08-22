# Contributing to ampliQC 🧬

Thank you for your interest in contributing to **ampliQC**! We welcome bug reports, feature requests, documentation improvements, and code contributions.

---

## Code of Conduct

Please be respectful and constructive in all issues, pull requests, and discussions.

---

## Development Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/LaBiOmicS/ampliQC.git
   cd ampliQC
   ```

2. Create a virtual environment and install development dependencies:
   ```bash
   python -m venv venv
   source venv/bin/activate
   pip install -e .[dev]
   ```

3. Run unit tests to verify your setup:
   ```bash
   pytest
   ```

---

## Submitting Pull Requests

1. Create a feature branch (`git checkout -b feature/my-feature`).
2. Make your changes and write unit tests in `tests/`.
3. Ensure all tests pass (`pytest`).
4. Commit your changes (`git commit -m "Add feature X"`).
5. Push to GitHub (`git push origin feature/my-feature`) and open a Pull Request.

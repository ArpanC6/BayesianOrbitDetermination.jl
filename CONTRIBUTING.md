# Contributing to BayesianOrbitDetermination.jl

We welcome contributions from astrodynamics researchers, flight dynamics engineers, and scientific machine learning developers!

## Development Workflow

1. **Fork and Clone**:
   ```bash
   git clone https://github.com/ArpanC6/BayesianOrbitDetermination.jl.git
   cd BayesianOrbitDetermination.jl
   ```

2. **Activate Environment**:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

3. **Run Test Suite**:
   ```julia
   using Pkg
   Pkg.test()
   ```

4. **Code Style**:
   - Follow standard SciML / Julia format guidelines.
   - Include docstrings for all exported functions.
   - Maintain unit test coverage for new dynamics or observation models.

5. **Submitting Pull Requests**:
   - Create a feature branch (`git checkout -b feature/my-feature`).
   - Ensure all tests pass.
   - Submit PR with clear description of physical models or algorithmic improvements.

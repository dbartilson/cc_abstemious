#  Development Guide

## Build and test

`cargo` is used for the build and test phases. If you install Rust via `rustup`, you should have `cargo` as well.

### Debug build

`cargo build` will build `cc_abstemious` in debug mode.

### Automated tests

```
cargo test
``` 

will build `cc_abstemious` in debug mode, then run the unit and integration tests.

`cargo test --release` will do the same, but in release mode.

## Documentation
```
cargo doc --no-deps
``` 

is run to generate the documentation. See the relevant GitHub Action which utilizes GitHub Pages. Note that the `index.html` in the `docs` directory is copied in to act as a redirect for that page.

## Dependency/SBOM Generation

[cargo cyclonedx](https://github.com/CycloneDX/cyclonedx-rust-cargo?tab=readme-ov-file) is used to generated the CyclondeDX JSON SBOM with this command:

```
cargo cyclonedx -f json
```

## GitHub Actions

There are two actions:

1. **On Merge Request** (to main branch): Run tests
2. **On Push** (to main branch): Run tests and deploy documentation to GitHub Pages.
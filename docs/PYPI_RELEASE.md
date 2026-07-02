# PyPI release checklist

This document describes the recommended release process for `structrace`.

## 1. Prepare the package metadata

Before every release, confirm:

- `pyproject.toml` has the intended `version`.
- `README.md` renders correctly on PyPI.
- `LICENSE` matches the intended project license.
- `python -m pytest` passes.

PyPI does not allow re-uploading the same version. If a release needs a fix,
increase the version, for example from `0.1.0` to `0.1.1`.

## 2. Build locally

From the repository root:

```bash
python -m pip install --upgrade build twine
Remove-Item dist -Recurse -Force
python -m build
python -m twine check dist/*
```

On macOS/Linux, replace the cleanup command with:

```bash
rm -rf dist
```

## 3. Test on TestPyPI

Create a TestPyPI account and API token, then upload:

```bash
python -m twine upload --repository testpypi dist/*
```

Install from TestPyPI in a clean environment:

```bash
python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ structrace
```

Smoke-test the package:

```bash
python -c "from structrace.watermark import embed_text, decode_text; print(embed_text, decode_text)"
python -m structrace --help
```

## 4. Publish to PyPI with Twine

After TestPyPI validation:

```bash
python -m twine upload dist/*
```

Users can then install the package with:

```bash
python -m pip install structrace
```

## 5. Publish with Trusted Publishing

The repository includes `.github/workflows/publish.yml` for PyPI Trusted
Publishing.

One-time PyPI setup:

1. Go to the PyPI project publishing settings.
2. Add a trusted publisher for GitHub.
3. Use repository `JLU-WangXu/Structrace`.
4. Use workflow file `.github/workflows/publish.yml`.
5. Use environment name `pypi`.

After setup, creating a GitHub Release or manually running the workflow can
publish the already-built distributions without storing a PyPI API token in
GitHub secrets.

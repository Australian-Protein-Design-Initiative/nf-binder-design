# Documentation

This site uses [MkDocs](https://www.mkdocs.org/) and the `readthedocs` theme, with [mike](https://github.com/jimporter/mike) for version management.

## Online Deployment

The online documentation is automatically deployed via a GitHub Action (`.github/workflows/docs.yml`). The action is triggered on pushes to `main`, `develop`, and version tags (either `v*` or unprefixed semver like `0.3.0`).

- Pushes to `main` deploy to the `/main/` directory and update the `latest` alias. `main` is also set as the default version.
- Pushes to `develop` deploy to the `/develop/` directory and update the `dev` alias.
- Pushes to version tags (e.g. `0.3.0`) are deployed as immutable versions (`/0.3.0/`).

## Local Development & Viewing

To view the documentation locally, complete with the version selector dropdown, you must deploy and serve the documentation using `mike` instead of a standard `mkdocs serve`.

### 1. Install prerequisites

Make sure you have the required packages installed:
```bash
pip install mkdocs mkdocs-material mike
```

### 2. Deploy your current branch locally

You need to use `mike deploy` to build your current branch into your local `gh-pages` branch before it can be served.

```bash
# E.g., if you are currently testing the 'develop' branch:
mike deploy --config-file docs/mkdocs.yml develop --update-aliases
```

If you want the root URL (e.g., `http://localhost:8000/`) to correctly redirect to a version, you'll need to fetch and deploy the default version, or set it localy:
```bash
mike set-default --config-file docs/mkdocs.yml latest
```

### 3. Start the Server

Start the local server using `mike serve` (which serves out of the `gh-pages` branch, meaning you will see the versioned dropdowns).

```bash
mike serve --config-file docs/mkdocs.yml
```

The documentation can then be viewed in your browser (e.g., `http://localhost:8000/`). Simply navigate to the specific version folder (e.g. `http://localhost:8000/develop/`) or the root to be redirected.

*(Alternatively, to quickly preview live-reloading changes on a single version while authoring without the versions dropdown, you can still use the standard `mkdocs serve --config-file docs/mkdocs.yml`).*

# how to use `mkdocs`

[`mkdocs`](https://www.mkdocs.org/) is a snazzy tool for making documentation easier. Zach learned about it from Christina Hedges' very useful tutorial [here](https://christinahedges.github.io/astronomy_workflow/notebooks/3.0-building/mkdocs.html). `mkdocs` is included as a `[develop]` dependency, so if you instaled this package via something like `pip install -e '.[develop]'`, you should already have it.

From within our repository directory....

I ran `mkdocs new .` which made a `docs/` directory and a `mkdocs.yml`

I edited the `docs/index.md` page to make it a more friendly introduction to the package.

I added notebooks in `docs/*.ipynb` explaining how to use some aspects of the code. Because mkdocs will automatically run the notebooks in the background, it's best to clear all output from the notebooks before saving/committing, so they don't take up unnecessary space in the repository.

I copied over and modified `mkdocs.yml` to define the navigation and basic settings of the documentation page.

I ran `mkdocs serve`, and woah, a live version of the docs appeared when I pointed my browser to http://127.0.0.1:8000/the-friendly-stars/. If I make and save changes to any of the documentation source, it reruns and soon updates what I see in the browser!

I ran `mkdocs gh-deploy` to deploy a pretty version of the docs up at `https://zkbt.github.io/the-friendly-stars/`. For the sake of not making the deployment `gh-pages` branch annoyingly large, add the `--no-history` option to erase the repository each time.

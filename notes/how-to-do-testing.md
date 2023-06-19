# how to use `pytest`

[`pytest`](https://docs.pytest.org/en/) is a snazzy tool for automatically testing that your code works. `pytest` is included as a `[develop]` dependency, so if you instaled this package via something like `pip install -e '.[develop]'`, you should already have it.

I include a folder called `tests/` inside my repository. Inside that folder I have modules called `test_*.py`, which contain functions that have the word `test` located somewhere in their name.

I periodically run `pytest` from the base repository directory. It searches through, finds any functions with `test` in their name, and runs them. If they work, it says so. If not, it doesn't!

Whenever I write code that want to be dependable, I add a test that I believe should run as long as the code is working. That way, as I changes things later on, I can always check that I haven't broken anything essential.

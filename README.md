# Code Examples
A place to share working examples of code useful to the lab.

# Structure of a Good Code Example
A good code example should have its own directory. 
Everything needed to run the code should be present in the directory (PDB files, python code) or clearly listed (conda packages).
There should also be a clear and detailed README.md file stating how the code can be used and outlining the theoretical or mathematical background to the code.
References for the theory, math, or code should be included in ACS format, including DOIs.
A user should be able to clone the repository and run the code that you've pushed within 15 minutes.

# Code Example Checklist
* Description of the goal of the code
* Description of how to use the code
* Description of the theoretical or mathematical background to the code
* References for the theory, math, or code with DOIs
* The code itself
* Files needed to run the code
* An example output

# Pushing and Pulling

1. Nativate to the repo on your local computer. 
2. `git pull`
3. Create a new branch using `git checkout -b [branch name]`.
4. Make the edits to your code examples.
5. `git status`
6. `git add` any files with changes or untracked files. 
7. `git commit -m "[descriptive explanation of what was done]"`
9. `git push --set-upstream origin [branch name]`
9. `git checkout main`  - switch back to the main branch.
10. Create a pull request. Go to [https://github.com/Daly-Lab-at-Haverford/code_examples](https://github.com/Daly-Lab-at-Haverford/code_examples). There should be a big green button at the top saying something like "Pull Request" Add a quick description of what you did then click create pull request.
11. Lab members will review your changes then approve your pull request, then delete your branch.
12. Delete your local copy of the branch with `git branch -d [branch name]`

# Cloning the repo to a new computer: 

1. `conda create -n gh_env`
   * You will be able to run git commands in any environment, but you will need to activate this `gh_env` to run gh commands.
2. `conda activate gh_env`
3. `conda install gh --channel conda-forge` 
4. `gh auth login`
    * Log into GitHub.com
    * HTTPS
    * Yes, authenticate Git with Github credentials
    * Use a web browser (don't forget to copy the one-time code)
5.  `git clone https://github.com/Daly-Lab-at-Haverford/code_examples.git`
6.  `conda deactivate`

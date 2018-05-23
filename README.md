# Analysis project template

This README contains detailed instructions on how to organize an analysis project. This repository also serves as a template repository you can clone as a starting point for a new analysis project.

## Components of a project

A data analysis project has four components:

* Raw data
* Metadata (descriptions of the data)
* Code
* Results (processed data, figures, *etc.*)

Which of these components belongs in a git repository? Only what we need under version control. Raw data and Results are usually too large, so they should not be version controlled and therefore should not be in the git repository. In contrast, metadata and code are typically small, text-based, and change frequently; these should therefore under version control in a git repository. Keeping the repository small improves efficiency, and then each person working on the project can clone it into a personal area because these duplicates are cheap and small. The larger components (Raw data and results) will instead reside in single copy on shared disk space, to be referenced by code or metadata in the git repository. 

Our project thus must be divided into two components: the version-controlled part and the non-controlled part, each with its own space and structure:

### 1. Version-controlled components: GitHub repository organization (metadata and code)

The two version-controlled components should each live under a subfolder. So, the simple structure would look like this:

```
github_repository
	/metadata
	/src
```

This is how we have organized this template repository. All code goes in `/src`. Descriptions, config files, and pointers to data go in `/metadata`. Pretty simple. The way to describe the project `metadata` is really important, too -- we want tools that can work with all our projects, so we have to define them in the same way. Therefore, we use a shared, standardized way of describing and configuring projects and sample annotation: [the PEP format](http://pepkit.github.io), which uses a `yaml` configuration file and a `csv` sample annotation sheet.

You should also think carefully about how you name this repository. See section below on [choosing a project name](#choosing-a-project-name).

### 2. Non-version-controlled components: Data organization (raw data and results)

Data and Results never get into the repository, so they must be stored elsewhere, in a shared disk area. They should be further divided into two areas due to the [the curse of enormity](http://databio.org/posts/curse_of_enormity.html) and different disk requirements for data and results:

1. Data is stored in a parent Data folder (either `$DATA` or `$MDATA`), subdivided into subfolders by data source (which are not required to have a one-to-one relationship with a project, as multiple projects may analyze the same data). The subfolder name therefore doesn't necessarily correspond to the project name, but instead, is more descriptive of the data itself.

2. Results are stored in the `$PROCESSED` folder, subdivided into a subfolder for each project. Because there *is* a one-to-one relationship between projects and result folders, **the results subfolder should be named the same as the git repository**, which should be a unique identifier for the project. Within that project subfolder, I like the concept of dividing Results into two parts: pipeline and analysis.  _Pipeline_ results are anything that takes you from raw, unaligned reads, to a processed data type, like RNA levels, methylation calls, or ChIP peaks. _Analysis_ is anything that happens afterwards. While the division is a bit ambiguous and probably not strictly necessary, I find it improves my conceptual approach to the project. Pipelines generally operate on an individual sample basis (or small subsets of samples, like 2-way comparison), and can be run as data is produced; they are generally relatively standardized by the community, with a series of steps primarily using public tools (plus some small internal supplements). Analysis is much more variable and less standardized, combining and comparing many samples to answer a particular question which varies by project. So a `$PROCESSED/project` folder would look something like this:
```
shared_data_folder
	/results_pipeline
		/sample1
		/sample2
	/analysis
		/fig1
		/fig2
		/tables
```


## How to start a new project:

(I have a script that automates this in [env](https://github.com/nsheff/env/blob/master/bin/newproject.sh))

1. Select a project name that fits with our naming style. See section below on [choosing a project name](#choosing-a-project-name).
2. Create a new repository on GitHub in our group organization (name `YOUR_PROJECT_NAME`)
3. Clone the `newproject` repository to establish folder structure (`git clone https://github.com/databio/newproject`)
4. Rename the local "newproject" folder to your new repo name (`mv newproject YOUR_PROJECT_NAME`).
5. Edit `metadata/project_config.yaml`. For consistency, the `output_dir` should have the same name as the project, so set `metadata.output_dir:` to `${PROCESSED}YOUR_PROJECT_NAME`. 
6. Edit this README.md file to briefly describe the project.
7. Delete the Makefile (this just lets me purge this history of this template so your new clone starts clean)
8. Update the remote of your local repo and push changes:

```
git remote set-url origin git@github.com:databio/YOUR_PROJECT_NAME.git
git push -u origin master
```

## Choosing a project name

To keep things simple and predictable, let's keep a correspondence between project name, GitHub repository name, and results folders. This means the project name should **match the output folder name exactly**. This 1:1 correspondondence makes it easier to keep track of things and find stuff when we're browsing around on disk, and also to build automated tools to look in results folders. Therefore, picking a name is pretty critical, because it's going to show up everywhere.

Let's also try to keep things consistent with our analysis project names. There's no right way, but there's the way we've been doing it, so we can just stick with that. Here are a few example project names:  `ews_patients`, `cphg_atac`, `microtest`, `ratbrain`.

Some style things:

* keep the project name all lowercase
* use underscores, not camelCase
* keep it relatively short, but long enough to capture the unique things about the project
* try to make it future-proof
* don't worry about it too much


## Project management guidelines for keeping GitHub repositories tidy

1. Put all source code in the `src/` folder.
2. Do not commit any results, data, figures, or binary files of any type into the repository.
3. Write [good commit messages](https://gist.github.com/nsheff/868f88bdca529e4a32377e8279dae826)
4. Keep the first line of commit messages [under 50 characters](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html) (use Details to elaborate)
5. Consider `--rebase` for git pull to [keep a cleaner history](http://git-scm.com/book/en/v2/Git-Branching-Rebasing)

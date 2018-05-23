clearhist:
	rm -rf .git
	git init
	git add .
	git commit -m "Initialize project template"
	git remote add origin git@github.com:databio/newproject.git
	git push -u --force origin master
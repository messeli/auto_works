# auto_works

The Auto demo files contain annotations, as well as the pdf's.  

# Extreme commenting guide:
The code mostly contains comments of the folowing format.

#| EXTREME COMMENTING LEGEND - self explained: 
#| 00  "# "        - a commented out line, either as a depot or as an option of code
#| 01  "#| "       - comment for the following lines
#| 02  "#|"        - comment at the end of a line
#| 03  "#|:"       - the above line/calculation is explained 
#| 04  "#| :"      - indentation of ":" suggests dt this is the explanation of "#| " (see 01 above)
#| 05  "#|::"      - explanation of an explanation of sort "#|:"; extend this w/ the idea in 04.
#| 06  "#|Nt:"     - a note not strictly related to the immediate preciding line, cld b deductions
#| 07  "#|XXX"     - a block of code starts here 
#| 08  "#|X_YYY"   - a block of code under the block "XXX"
#| 09  "#|XY_ZZZ"  - a block of code under the block "X_YYY", d reference to d X block cn be omitd
#| 10  "#|XXX END" - a block of code ends here
#| 11  "#|..."     - continued line from either "#|:" or "#|::" or "#|XXX" or "#|Nt:" or "#|..." 
#| 12  "#| ..."    - continued line to the 1 above (also see 04 above for another case) 
#| 13  "xxxYyyZzz" - a term generated on the spot, could save space; could be a variable name, too
#| 14  shortHands  - d:the | s:is | t:to | b:be | r:are | f:of | u:you 
#| 15                ds:this | dt:that | dn:then or than | cn:can | cld:could | fr:for 
#| 16                dss:this is | dts:that is   
#| 17                recomd:recommended | verbd: pastParticipleForm of a verb 
#| 18                J or G: -ing
#| 19                w/:with | w/o:without 
#| 20                others...
#|
#|EEE Examples
#| See some examples 
#| ...on 09.
#| :As dss a self-explained legend
some code = some more code 
#|:dss some code here. 
#| ::Self-explained here means it has its examples in it.
#|EEE END Examples
#|==================================================================================================


Git repo start guide:
1)	Setup Git 
a.	Ubuntu: sudo apt-get git 
b.	Windows: https://git-scm.com/download/win
2)	git config --global user.name ‘Your Name’
3)	git config --global user.email ‘youremail@domain’
4)	git config –global credential.helper store
a.	This will keep you signed in for GitHub when pulling or pushing.
5)	Navigate  to a directory to store your local Git repository. Then, 
a.	Ubuntu: Right click > Open Terminal
b.	Windows : Right click > Git Bash Here 
6)	git init
7)	git branch -M main
a.	Default branch name was “master”, this will rename it to “main” to conform to the GitHub repo branch name. 
b.	If it is not set to main, when pushing it will generate a new remote repo branch with local branch’s name.
8)	git pull https://github.com/messeli/auto_works.git
9)	Enter your GitHub username and password (Maybe I must specify your Git username as collaborators)
10)	git remote add origin https://github.com/messeli/auto_works.git 
a.	Remote repo to push to and pull from  is named origin by default.
11)	Do changes, add files etc …
12)	git add .
13)	git commit -m ‘Your Comment.’
14)	git push -u origin main
15)	After this first push you can simply run git push to record your commits in the remote repo.
16)	Remember always to do git pull before working. 
a.	Someone might have changed something overnight. 

Additional notes: 
1)	Git restore <file> 
a.	This will discard your unstaged changes. It will bring to the working directory what it was in the staging area. If you have not staged anything, then it will delete it. (!?)
2)	Git checkout <file>
a.	Will lose your working directory content similar to the way above.
3)	Some other undoing methods might be destructive.
4)	Technically anything that is committed should be somehow accessible. 
a.	Though you can still lose them if you do some commit modifications wrong. 
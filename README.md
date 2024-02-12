# auto_works

The Auto demo files contain annotations, as well as the pdf's.  

# Extreme commenting guide:
The code mostly contains comments of the folowing format.

#| EXTREME COMMENTING LEGEND - self explained: 
#| 00  "# "        - commentD line, as depot or an option of code
#| 01  "#| "       - comment describJ d followJ lines
#| 02  "#|"        - comment at d end f a line
#| 03  "#|:"       - comment fr d above line
#| 04  "#| :"      - indentaTn f ":" says dt dss d explanaTn f "#| "
#| 05  "#|::"      - explanaTn f an explanaTn f sort "#|:"
#| 06  "#|Nt:"     - note nt strictly relatD t precidJ line; deducTn
#| 07  "#|:Nt:"    - note strictly relatD t precidJ line: cd b deducTn
#| 08  "#|XXX>"    - block f code starts, also cl b "#|___XXX>"
#| 09  "#|X_YYY>"  - block f code under d block "XXX" 
#| 10  "#|XY_ZZZ>" - block f code under d block "X_YYY"
#| 11  "#|XXX."    - block f code ends here, old versions:"XXX END"
#| 12  "#|..."     - continD frm ["#|:","#|::","#|Nt:","#|..."] or other
#| 13  "xxxYyyZzz" - term generatD onDspot, cd save space; cd b a var
#| 14  shortHands  - d:the | s:is | t:to | b:be | r:are | f:of | u:you 
#|                   ur:your | dnt: don't | dsnt: doesn't
#|                   ds:this | dt:that | dn:then or than
#|                   cn:can | cd:could | fr:for 
#|                   dss:this is | dts:that is   
#|                   verbD: verb3 | verbJ,G: -ing | Tn,Sn: -tion,-sion 
#|                   w/:with | w/o:without 
#|                   others...
#|
#|___EEE> Examples
#| See some examples 
#|...on 09.
#| :As dss a self-explainD legend
#| ::Self-explainD here means it has its examples in it.
#|___EEE. Examples
#|======================================================================


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

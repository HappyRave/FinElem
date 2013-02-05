#
#  Makefile for MAC OS (cocoa)
#
#     make all   (construction de l'executable)
#     make clean (effacement des fichiers objets et de l'executable)
#
#  A adapter en fonction des ordinateurs/environnements 
#  Compilateur, edition de liens, 
#
GLFWINCLUDE  = /Users/me/glfw277/include
GLFWLIBRARY  = /Users/me/glfw277/lib/libglfw.a
#
CC       = gcc  
LD       = gcc
CFLAGS   = -O3 -Wall
LFLAGS   = -Wall -O3
LIBS     = -lm 
#
PROG     = myFem
LISTEOBJ = \
  homework.o  main.o
# ATTENTION... aucun caractere apres le caractere de continuation "\"
#
# compilation
#
.c.o :
	$(CC) -c  $(CFLAGS) -o $@ $<
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# dependances
#
all        : $(PROG)
homework.o : homework.c 
main.o     : main.c 
#
# edition de lien
#
$(PROG) : $(LISTEOBJ)
	$(LD) -o $(PROG) $(LFLAGS) $(LISTEOBJ) $(LIBS)
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# effacement des fichiers intermediaires
#
clean :
	rm -vf $(PROG) $(LISTEOBJ) core a.out
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# ATTENTION... il faut une ligne vide a la fin du fichier.



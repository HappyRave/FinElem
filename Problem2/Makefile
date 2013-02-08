#
#  Makefile for MAC OS (cocoa)
#
#     make all   (construction de l'executable)
#     make clean (effacement des fichiers objets et de l'executable)
#
#  A adapter en fonction des ordinateurs/environnements 
#  Compilateur, edition de liens, 
#
GLFWINCLUDE  = /usr/local/include/GL/include
GLFWLIBRARY  = /usr/local/include/GL/libglfw.a
#
CC       = gcc  
LD       = gcc
CFLAGS   = -I$(GLFWINCLUDE) -O3 -Dgraphic -Wall -O0 -g -Qc--enable-checking
LFLAGS   = -Wall -O3 -framework Cocoa -framework OpenGL -framework IOKit 
LIBS     = $(GLFWLIBRARY) -lm 
#
PROG     = myFem
LISTEOBJ = \
  homework.o   fem.o   glfem.o  main.o
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
homework.o : homework.c fem.h
fem.o      : fem.c fem.h
glfem.o    : glfem.c fem.h glfem.h
main.o     : main.c fem.h glfem.h
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



#It only makes sense to create all these CFLAGS and LDLIBS if we intend to use them and they weren't
FLAGS = -o -w -std=gnu99 
#fo the gsl, it must be written as $(shell location/gsl-config --"CFALGS/LDLIBSpart") just written gsl-config is unfortunately not enough 
CFLAGS += 'gsl-config --cflags' # you cant use ' ' use $() instead as in the next line 
LDLIBS += $(gsl-config --libs)

# The ";" should not be there for any of them

all: output ; cat $@		# additional commands are to be placed in the next line  also the $@should not be there
	output: main : $@ > $<		#We are printing out our result so we are creating a txt file. We must say this by written output.txt. Again additional commands must be place in the next line after and indent not this weirdness with :. $@ and $< have been placed the wrong place. $@ is the target "output" $< is main not the over way around. You also need to start the command with ./ 

#the creation of the main is missing and since we are using headers like gsl we need to write the command explicite and put in (CFLAGS) and (LDLIBS) in it


clean: ; $RM main output	#again commands are placed in the next line, and $RM should be written as $(RM)



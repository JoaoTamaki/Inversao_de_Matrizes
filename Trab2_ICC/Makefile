    CC     = gcc -g
    CFLAGS =
    LFLAGS = -lm -O3 -mavx -march=native

      PROG = invmat
      OBJS = utils.o \
             sislin.o \
             Metodos.o \
             $(PROG).o

.PHONY: limpa faxina clean distclean purge all

%.o: %.c %.h utils.h sislin.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak
	@rm -f *.o core a.out
	@rm -f $(PROG)

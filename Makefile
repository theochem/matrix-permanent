CC ?= gcc

CFLAGS = -fPIC -Wall -Wextra -O2 -g

SRCS = permanent/permanent.c
OBJS = $(SRCS:.c=.o)

.PHONY: all
all: libpermanent.so

libpermanent.so: $(OBJS)
	$(CC) -shared -o $@ $^

$(SRCS:.c=.d):%.d:%.c
	$(CC) $(CFLAGS) -MM $< >$@

include $(SRCS:.c=.d)

.PHONY: clean
clean:
	rm -f ${TARGET_LIB} ${OBJS} $(SRCS:.c=.d)

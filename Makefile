CC = csc
SRC = src
BIN = bin
TEST = test
TEST_FILES = test-files
TEST_OUTF = $(TEST_FILES)/rec.bam.bin-cov.samtools-depth-aa.graph-lines.cov-plot.pdf

MAIN = bin-cov

SAMTOOLS = `which samtools`

# TEST_SRC := $(wildcard $(TEST)/*.scm)
# TEST_OBJS := $(TEST_SRC:.scm=.o)
# TEST_BINS := $(TEST)/*.bin

# SRC_FILES := $(wildcard $(SRC)/*.scm)
# SRC_OBJS := $(SRC_FILES:.scm=.o)

.PHONY: all
.PHONY: clean
.PHONY: install
.PHONY: test
.PHONY: test_clean

all: $(MAIN)

clean:
	-rm -r $(SRC)/*.o ./$(MAIN) $(TEST_OUTF) $(TEST_FILES)/*.bin-cov.*

install:
	cp ./$(MAIN) /usr/local/bin

test: $(MAIN)
	bin/bin-cov $(SAMTOOLS) $(TEST_FILES)/rec.bam $(TEST_FILES)/names.txt

# test: compile_tests
# 	for test in $(TEST_BINS); do \
# 	  $$test; \
# 	done

# compile_tests: $(SRC_OBJS) $(TEST_OBJS)
# 	$(CC) -o $(TEST)/pfa-test.bin $(SRC)/pfa.o $(TEST)/pfa-test.o

# %.o: %.scm
# 	$(CC) -c $?

$(MAIN):
	mkdir -p bin
	$(CC) -c $(SRC)/bin-cov.scm
	$(CC) -o $(BIN)/$@ $(SRC)/bin-cov.o

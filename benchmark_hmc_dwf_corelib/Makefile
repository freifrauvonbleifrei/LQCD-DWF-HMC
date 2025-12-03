.SUFFIXES: .cpp .h .o

BRIDGE_INSTALL_PATH = ../build

include $(BRIDGE_INSTALL_PATH)/make.inc

headers =
sources = hmc_dwf_eo.cpp
objects = $(sources:%.cpp=%.o)

program = hmc_dwf_eo.elf

all: msg $(program)

$(program): $(objects)
	$(LD) $(LDFLAGS) -o $@ $(objects) $(LDLIBS)

%.o: %.cpp $(headers)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

msg:
	@cat $(BRIDGE_INSTALL_PATH)/message.txt

clean:
	rm -f $(program) $(objects)


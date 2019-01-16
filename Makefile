compiler = g++
flag = -fopenmp
excutable = PandasCommute.x
headfile = $(shell find code/ -name "*.h")
src = $(shell find code/ -name "*.cpp")
obj = $(src:%.cpp=%.o)

$(excutable) : $(obj)
	g++ -o $(excutable) $(obj) $(flag)

%.o : %.cpp $(headfile)
	$(compiler) -c $< -o $@ 

clean:
	rm -rf $(obj) $(excutable)

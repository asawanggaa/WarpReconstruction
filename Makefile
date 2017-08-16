CPPFLAGS = -I/usr/include -I/usr/local/include
LDLIBS = -lglut -lGLU -lGL -lGLEW -lm \
	-lopencv_calib3d -lopencv_contrib -lopencv_core -lopencv_features2d \
	-lopencv_flann -lopencv_gpu -lopencv_highgui -lopencv_imgproc -lopencv_legacy \
	-lopencv_ml -lopencv_objdetect -lopencv_ocl -lopencv_photo \
	-lopencv_stitching -lopencv_superres -lopencv_ts -lopencv_video -lopencv_videostab
LDFLAGS = -std=c++11 -O2

# CPPFLAGS = 'pkg-config --cflags opencv'
# LDLIBS   = 'pkg-config --libs opencv'

objects = TestMain.o ShapeWarp.o PredictedPath.o ConfigurableKCurve.o

TestMain : $(objects) 
	clang++ $(LDFLAGS) -o TestMain $(objects) $(CPPFLAGS) $(LDLIBS)
TestMain.o : TestMain.cpp PredictedPath.h ShapeWarp.cpp ShapeWarp.h TestFunction.hpp 
	clang++ $(LDFLAGS) -c TestMain.cpp $(CPPFLAGS)
ShapeWarp.o : ShapeWarp.h PredictedPath.h ShapeWarp.cpp
	clang++ $(LDFLAGS) -c ShapeWarp.cpp $(CPPFLAGS)
PredictedPath.o : PredictedPath.h ConfigurablePath.h PredictedPath.cpp
	clang++ $(LDFLAGS) -c PredictedPath.cpp $(CPPFLAGS)
ConfigurableKCurve.o : ConfigurablePath.h ConfigurableKCurve.cpp
	clang++ $(LDFLAGS) -c ConfigurableKCurve.cpp $(CPPFLAGS)

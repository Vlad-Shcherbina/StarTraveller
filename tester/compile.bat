mkdir build
javac -d ./build StarTravellerVis.java
cd build
jar -cvfe tester.jar StarTravellerVis *
cd ..
move build\tester.jar .

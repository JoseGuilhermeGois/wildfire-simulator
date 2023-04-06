# IMfire

In order to simulate the fire spread, the terrain we're evaluating will be treated as a 2d matrix of cells. For easier tranlation to the real world certain analogies will be applied for naming conventions.

* cell >> element
* matrix >> topography plane


## Description of some parts of the code:
### Config Folder

Configuration folder is to read the files needed for the simulation.
* Environment
* Fuel
* Landscape


1. BusinessConfigProcessor(business_processor.py) 
    - Abstract class that read the file path for every file type.
    - This class must not be instantiated.
    - All subclasses have to implement process.
2. landscape folder 
   1. landscape.py
       - (Landscape) This is where we save the values that we get from reading the data inside landscape folders.
       - We have to create two different classes for the (Shape) and for the (Location) to then save these values on the (Landscape) class.
   2. landscape_processor.py
       - This is where we read the file.
       - First we read the shape, then location and then element size.
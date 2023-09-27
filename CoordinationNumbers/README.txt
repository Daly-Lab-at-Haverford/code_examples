This will folder will contain a code example of calculating coordination numbers. A coordination number calculation will give you the specific number of one atom around another atom. 

To provide a better picture, I will be modeling a coordination number calculation of finding the amount of water molecules around another water with oxygen being the reference atom. The reference atom is just where in your molecule do you want the calculation to begin since coordination numbers calulate radially beginning at that "reference atom". 

To also make it easier for you to navigate through the files, I have included everything you will need to observe how to understand the code, or to replicate it within your own research. 

Please navigate to the 'Materials' folder within the 'CoordinationNumbers' folder to get a glimpse of what you will need to. 

I have also changed the directory of these files from my own work into a folder under 'kndiaye' titled 'code_materials' and then 'coordination_numbers' to see the exact files I'm working with. When working with coordination numbers, you want to have:
    - a trajectory to analyze (ex. output.h5)
    - an RDF calculation of the atoms you want to get a read on
        - if you do not know how to calculate a RDF, please go to the 'RDF' folder in the 'code_examples' repo
    - an output.dat file of your trajectory results 
        - you are not going to add the output.dat file into your code, BUT it is good to have on hand because there are specific values that you are going to need from it.
    
After reviewing the 'Materials' folder, please navigate to the 'Code' jupyternotebook. Here is where the code is and there are comments in each of the cells for better understanding.

(If you accidentally edit my example materials in the 'code_materials', that is completely ok! They are copies of my original work which are safely in separate folders. Just let me know so I can restore them to their original state)

The directory of the 'code_materials' to use as reference as you go through the tutorial is:

`/homes/kndiaye/code_materials/coordination_numbers`
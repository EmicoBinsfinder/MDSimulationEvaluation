'''
File to edit dihedrals for LOPLS
Do so by changing format of OPLS settings file
to allow use with hybrid dihedral style
'''

import os

def modify_dihedral_coeff_lines(input_file, output_file):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if 'dihedral_coeff' in line:
                    # Splitting the line by spaces to insert "multi/harmonic" after "dihedral_coeff"
                    parts = line.split()
                    # Modify the line if it starts with 'dihedral_coeff'
                    parts.insert(2, 'opls')
                    # Join the parts back into a single line and write to the output file
                    modified_line = ' '.join(parts)
                    print(line)
                    print(modified_line)
                    outfile.write(modified_line + '\n')
                else:
                    # Write the unmodified line to the output file
                    outfile.write(line)
    except FileNotFoundError:
        print(f"Error: The file {input_file} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
input_filename = 'system.in.settings'  # Replace with your actual input file name
output_filename = 'output.txt'  # Desired output file name

modify_dihedral_coeff_lines(input_filename, output_filename)

        


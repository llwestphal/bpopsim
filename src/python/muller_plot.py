# Code by Aaron Reba
# 12-09-24

import argparse
import tempfile
import os

import matplotlib.pyplot as plt
import Image
import numpy

def check_adjacent(text_data, x, y):
    checking = []
    if x != 0:
        checking.append((x - 1, y))
    if x != text_data.shape[0] - 1:
        checking.append((x + 1, y))
    if y != 0:
        checking.append((x, y - 1))
    if y != text_data.shape[1] - 1:
        checking.append((x, y + 1))
    
    adjacent = []
    
    for check in checking:
        if text_data[x][y] != text_data[check[0]][check[1]]:
            adjacent.append(check)
    return adjacent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='file_name')
    parser.add_argument('-o', '--output', action='store', dest='save_name')
    parser.add_argument('-x', '--width', action='store', dest='width_expansion', default=1, type=int)
    parser.add_argument('-y', '--height', action='store', dest='height_expansion', default=1, type=int)
    results = parser.parse_args()
    file_name = results.file_name
    save_name = results.save_name
    width_expansion = results.width_expansion
    height_expansion = results.height_expansion
    
    #make colors
    color_list = []
    for r in xrange(0, 256, 127):
        for b in xrange(0, 256, 127):
            for g in xrange(0, 256, 127):
                color_list.append((r, g, b))
    color_list.remove((0, 0, 0))
    color_list.remove((127, 127, 127))
    color_list.remove((254, 254, 254))
    
    #temp = color_list[:3]
    #color_list = temp
    
    text_data = numpy.genfromtxt(file_name, dtype=None)
    
    width, height = text_data.shape
    
    plot = Image.new('RGB', (height * width_expansion, width * height_expansion))
    
    unique_types = []
    unique_locations = {} #indexed as unique
    #yields a list of 3 lists:
    #list of coords of that type
    #list of surrounding coords that are not that type
    #list of all coords belonging to this type
    
    #find unique types
    for x in xrange(width):
        for y in xrange(height):
            this_type = text_data[x][y]
            if this_type not in unique_types:
                unique_types.append(this_type)
                unique_locations[this_type] = [[], set(), []]
            adjacent = check_adjacent(text_data, x, y)
            if adjacent:
                unique_locations[this_type][0].append((x, y))
                unique_locations[this_type][1].update(adjacent)
            unique_locations[this_type][2].append((x, y))
    
    #while all types aren't touching
    retry = True
    while retry:
        color_map = {} #indexed as (x, y), yields rbg
        color_index = 0
        retry = False
        
        #make color map
        for unique in unique_types:
            
            #if color_index loops back around to start_index, all colors
            #have been checked and the current configuration is impossible.
            start_index = color_index
            #looping colors
            for i in xrange(len(color_list)):
                color = color_list[color_index]
                #check to see if this type has any neighbors that are the
                #same color
                neighbor_position_list = unique_locations[unique][1]
                for neighbor_position in neighbor_position_list:
                    if neighbor_position in color_map:
                        if color_map[neighbor_position] == color:
                            #neighbor with same color found,
                            break
                else:
                    #a clean loop exit.
                    #no neighboring non-identical positions found that share
                    #same color. this is a good color choice. go on to the next
                    #set of uniques
                    
                    #load colors into color_map
                    type_position_list = unique_locations[unique][2]
                    for type_position in type_position_list:
                        color_map[type_position] = color
                    break
                
                color_index += 1
                if color_index == len(color_list):
                    color_index = 0
                
                if color_index == start_index:
                    #all colors have been tried, this combination will not
                    #work
                    retry = True
                    print 'Need more colors.'
                    return
                    break
            else:
                #clean exit. update color.
                color_index += 1
                if color_index == len(color_list):
                    color_index = 0
                
            if retry:
                break
            
        else:
            #a clean loop exit.
            #a good combination has been found.
            pass
    
    #draw image
    for y in xrange(height):
        for x in xrange(width):
            for h in xrange(height_expansion):
                for w in xrange(width_expansion):
                    plot.putpixel(((height - y - 1) * width_expansion + w, x * height_expansion + h), color_map[(x, y)])
    
    plot = plot.rotate(90)
    
    plot.save(save_name)

main()
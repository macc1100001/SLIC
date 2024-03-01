#!/bin/bash

valgrind -s --show-leak-kinds=all --leak-check=full ./slic lena_std.tif

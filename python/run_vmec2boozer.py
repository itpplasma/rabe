if __name__ == "__main__":
    import sys

    from libneo import BoozerFile

    input_file = sys.argv[1]
    output_file = input_file + ".bc"

    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        output_file = "field.bc"

    if len(sys.argv) >= 4:
        uv_grid_multiplier = int(sys.argv[3])
    else:
        uv_grid_multiplier = 6

    b = BoozerFile(filename="")
    b.convert_vmec_to_boozer(
        filename=input_file, uv_grid_multiplicator=uv_grid_multiplier
    )
    b.write(output_file)

search_dir='/Users/jgoldman/Desktop/chp3/ltr_nbr/'
for f in "$search_dir"/*
do

    x="$(basename "$f" .tif)"
    b="$(basename "$f")"
    echo "$b"
    echo "$x"
    n=""$x".nc"
    
    gdal_translate -of NetCDF $b $n

done


dir=$1
i=0
for filename in "$dir"\/*.pickle; do
        echo "$filename"
	cp "$filename" "$dir"/"$i".in
        i=$((i+1))
done

for year in $(seq 1992 2014)

do
	echo 'processing '$year

	echo 'downloading file...'
	if [ -f "sanger_$year.html" ]; then
		echo "sanger_$year.html already downloaded."
	else
		curl -s https://www.sanger.ac.uk/research/publications/$year.html > sanger_$year.html
	fi

	echo 'cleaning and parsing file...'
	sed -n -e "47,$(($(wc -l < sanger_$year.html) - 80))p" sanger_$year.html | grep -B 3 -E "publication.*>Cell|publication.*>Nature|publication.*>Science" | grep -v reviews | grep -v Georgetown | grep "author" | sed 's_<span class="author">__g' | sed 's_</span>__g'| sed 's_<p class="authors">__g' | sed 's_</p>__g'| sed -e 's/, /\
/g' | sed -e 's/ and /\
/g' | awk 'NF<=2' | grep -v "&" | sort | sed 's/^ *//' | uniq -c | awk '$1 >= 1 {print $2$3"\t"$1}' > big3_$year.txt 

	echo 'cleaning and parsing file...'
	sed -n -e "47,$(($(wc -l < sanger_$year.html) - 80))p" sanger_$year.html | grep -B 3 -E "publication" | grep -v reviews | grep -v Georgetown | grep "author" | sed 's_<span class="author">__g' | sed 's_</span>__g'| sed 's_<p class="authors">__g' | sed 's_</p>__g'| sed -e 's/, /\
/g' | sed -e 's/ and /\
/g' | awk 'NF<=2' | grep -v "&" | sort | sed 's/^ *//' | uniq -c | awk '$1 >= 1 {print $2$3"\t"$1}' > all_publications_$year.txt 

	### explaination of the code

	# -> remove first 47 and last 80 lines
	#sed -n -e "47,$(($(wc -l < sanger_$year.html) - 80))p" sanger_$year.html

	# -> Extract lines with Cell, Nature or Science in the journal name, and grab the 3 lines before it
	#grep -B 3 -E "publication.*>Cell|publication.*>Nature|publication.*>Science" sanger_parsed.txt 

	# -> remove reviews, and remove some polluting cell journal name
	#grep -v reviews | grep -v Georgetown 

	# -> grab only the line containing the list of author names
	#grep "author" 

	# -> remove all the HTML formatting
	#sed 's_<span class="author">__g' | sed 's_</span>__g'| sed 's_<p class="authors">__g' | sed 's_</p>__g'

	# -> replace ", " and " and " characters with a new line (literally escapes the new line to keep sed running)
	#sed -e 's/, /\
	#/g' | sed -e 's/ and /\
	#/g'

	# -> throw out sporadic sentences that somehow appeared
	# awk 'NF<=2'
done

echo 'extracting faculty members'
if [ -f "faculty.html" ]; then
		echo "faculty.html already downloaded."
	else
		curl https://www.sanger.ac.uk/research/faculty/ > faculty.html 
	fi

grep h3 faculty.html | grep -v div | grep -v links | sed 's_</h3>__g'| sed 's_<h3>__g'| sed 's_-.*__g' | sed 's/^ *//' | awk '{print $2}' | grep -v Deloukas | grep -v Carter | grep -v Futreal | grep -v Peltonen > faculty.txt


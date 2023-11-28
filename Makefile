all:
	g++ -O3 -Wall --std=c++17 faiss2simple.cpp -o faiss2simple -ltbb
	gcc -O3 -Wall decoder.c -o decoder -lm 
	gcc -O3 -Wall encoder.c -o encoder -lm  
	gcc -O3 -Wall quantize.c -o quantize -lm

clean:
	rm faiss2simple
	rm decoder
	rm encoder
	rm quantize

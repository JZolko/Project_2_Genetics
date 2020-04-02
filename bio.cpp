#ifndef BIO_H
#define BIO_H

#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
using std::cout; using std::endl;
#include <vector>


/*
This function should return true if and only if
every character in the input is one of ATCG.
*/
bool IsValidDNASequence(const std::string & input){

  for(auto i = 0; i < static_cast<int>(input.size()); i++){

    if( input[i] == 'C' || input[i] == 'A' or input[i] == 'T' or input[i] == 'G'){}

    // If any of the letters arent good, return false
    else{
      return false;
    }
  }
  return true;
}


/*
This function should calculate the reverse complement DNA sequence.

The first argument is the sequence, the second argument is a pointer to
an empty string, which you should modify to store the result.

This is obtained by reversing the input sequence and swaping each
nucleotide/letter with it's complement:
A <-> T
C <-> G

Example:
input = AAATTCGGGG
reverse = GGGGCTTAAA
reverse complement = CCCCGAATTT
*/

void GetReverseComplementSequence(const std::string & input,  std::string * output){
  std::string normal = "ATCG";
  std::string normal_comp = "TAGC";
  size_t found = 0;

  //goes over it backwards to save from reversing the actual string
  for(int i = static_cast<int>(input.size()) - 1; i >= 0; i--){


    //Gets the index of a char in the string and adds the
    //compliment from the other string
    found = normal.find(input[i]);
    *output += normal_comp[found]; //adds opposite into string


  }
}
/*
This function should return the RNA transcript from a DNA sequence.

A RNA transcript is the reverse complement of the DNA sequence, but RNA
has U (uracil) instead of T (thiamine).

Make sure you don't have any redundant code.v.push_back("");
*/

std::string GetRNATranscript(const std::string & input){
  std::string output;
  std::string normal = "AUTGC";
  std::string normal_comp = "UAACG";
  size_t found = 0;

  //Iterates through string backwards to save from reversing string
  for(int i = static_cast<int>(input.size()) - 1; i >= 0; i--){

    //Gets the index of the char in the normal string
    //then adds the complement to the output string
    found = normal.find(input[i]);
    output += normal_comp[found];

  }
  return output;
}


/*
This function should return a vector of vector of strings with each possible RNA
reading frame from the given DNA sequence.

There are three possible reading frames (because the genetic code has three
nucleotides per amino acid) in each direction (you can also transcribe DNA in
the reverse complement direction, called the antiparallel strand).

Order the sequences like so:
1: Original (0 offset)
2: Original (1 offset)
3: Original (2 offset)
4: Antiparallel (0 offset)
5: Antiparallel (1 offset)
6: Antiparallel (2 offset)

With in the input sequence of: AATTCCCGAAA
Original RNA transcript = UUUGCCCAAUU
Antiparallel RNA transcript = AAUUCCCGAAA

The offsets (starting at pos 0, 1, and 2) of the two RNA transcripts
UUUGCCCAAUU
UUGCCCAAUU
UGCCCAAUU
AAUUCCCGAAA
AUUCCCGAAA
UUCCCGAAA

Instead of returning a vector of 6 strings, break each string into a vector
of length 3 strings (called codons) These codons will be useful
for the next translation step.

UUUGCCCAAUU -> {"UUU", "GCC", "CAA"}
// drop any remaining letters that don't fill a codon
UUGCCCAAUU -> {"UUG", "CCC", "AAU"}
UGCCCAAUU -> {"UGC", "CCA", "AUU"}
AAUUCCCGAAA -> {"AAU", "UCC", "CGA"}
AUUCCCGAAA -> {"AUU", "CCC", "GAA"}
UUCCCGAAA -> {"UUC", "CCG", "AAA"}
*/

std::vector<std::vector<std::string>> GetReadingFramesAsCodons(const std::string & input){

  std::vector<std::vector<std::string>> output;
  std::vector<std::string> transcriptCodons;
  std::vector<std::string> v;
  std::string NewInput;

  NewInput = GetRNATranscript(input);

  std::string NewNewInput = NewInput;

  int offset = 0;
  int cycle = 0;
  //I feel like i could've bnoken this next
  //atrocity into 2 different loops instead of doing something thats O(n^3)

  //loops twice to factor in regular and antiparallel
  while(cycle < 2){

    // accounts for offsets in indexing
    while(offset < 3){
      v.push_back("");

      //cycles trhough based on offset
      for(int i = offset; i < static_cast<int>(NewInput.size()); i++){

        // adds codons of size 3 to a string vector
        if(v.back().size() < 3){
          v.back() += NewInput[i];
        }
        else{
          v.push_back("");
          v.back() += NewInput[i];
        }
      }

      offset += 1;

      // adds only codons of size 3 into string vector
      for (int i = 0; i < static_cast<int>(v.size()); i++) {
        if(v[i].size() == 3){
          transcriptCodons.push_back(v[i]);
        }
      }

      // adds the string vector of codons into the
      // vector of string vectors
      output.push_back(transcriptCodons);

      //clears to avoid any mixmatching of codons
      transcriptCodons.clear();
      v.clear();

    }

    cycle += 1;

    // gets the antiparallel to turn into codons
    NewInput = GetRNATranscript(NewNewInput);

    // resets the offset for the new transcript
    offset = 0;
  }
  return output;
}



/*
This function translates/converts a vector<string> (vector of codons) into a
string of amino acids using the genetic code
(see https://en.wikipedia.org/wiki/Genetic_code).

For example, the codons:
{"UUU", "GCC", "CAA"}
translates to:
F (Phenylalanine), A (Alanine), Q (Glutamine)
abreviated:
FAQ

http://www.soc-bdr.org/rds/authors/unit_tables_conversions_and_genetic_dictionaries/genetic_code_tables/


To make your lives easier, here's a list of the possible codons:
"GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",
"AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
"GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
"UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
"CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
"ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
"GUG", "UAG", "UGA", "UAA"

And there corresponding amino acids ("*" represents STOP codons,
more on them later):

"A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D",
"C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
"I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
"P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
"Y", "V", "V", "V", "V", "*", "*", "*"
*/
std::string Translate(const std::vector<std::string> & codon_sequence){

  std::vector<std::string> codons = {"GCU", "GCC", "GCA", "GCG", "CGU",
  "CGC", "CGA", "CGG", "AGA", "AGG",
  "AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
  "GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
  "UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
  "CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
  "ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
  "GUG", "UAG", "UGA", "UAA"};

  std::vector<std::string> aminos = {"A", "A", "A", "A", "R", "R", "R",
  "R", "R", "R", "N", "N", "D", "D",
  "C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
  "I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
  "P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
  "Y", "V", "V", "V", "V", "*", "*", "*"};

  std::string out;

  //cout << "-------------------------"<< endl;
  //for(int i = 0; i < aminos.size(); i++){
  //  cout << aminos[i] << "  |  " << codons[i] << endl;
  //}

  //https://thispointer.com/c-how-to-find-an-element-in-vector-and-get-its-index/
  //auto found = 0;

  for (size_t i = 0; i < (codon_sequence.size()); i++) {

      //gets a pointer to a codon in the codon vector
      auto found = std::find( codons.begin(), codons.end(), codon_sequence[i] );

      // the "if" is probably redundant considering
      // all combos exist, but whatever
      if(found != codons.end()){

        //gets the index of the codon in the vector
        auto place = std::distance(codons.begin(), found);

        // adds the respective amino acid into the output/
        // string based on index
        out += aminos[place];
      }

  }

  return out;
}

/*
This function takes a DNA sequence and returns the longest possible
amino acid sequence / protein that is encoded by that sequence
(open reading frame). A valid open reading frame begins with thev.push_back("");
codon AUG (the amino acid, Methionine (M)) and runs until a stop codon (*)
is encountered. There may be multiple open reading frames in a sequence, and
you need to check all six reading frames in order given by
get_reading_frames_as_codons. If there are ties for longest, favor the first
one found.

Return the longest open reading frame as an amino acid sequence. It must start
with an 'M' and end with a '*' with no other '*''s within.
*/
std::string GetLongestOpenReadingFrame(const std::string & DNA_sequence){

  std::vector<std::string> v;
  std::string translated;
  std::string out;

  const std::vector<std::vector<std::string>> rfac = GetReadingFramesAsCodons(DNA_sequence);

  //to account for each list of codons
  for(int i = 0; i < static_cast<int>(rfac.size()); i++){
    translated = Translate(rfac[i]);


    v.push_back("");

    // goes through each char
    for (int i = 0; i < static_cast<int>(translated.size()); i++) {


      // starts adding to string vector if it goes over an M
      if(translated[i] == 'M'){

        if( v.size() > 0 && v.back().size() > 0){
          v.back() += translated[i];
        }
        else{

          v.push_back("");
          v.back() += translated[i];
        }
      }

      // stops adding and makes an empty string at
      // the back of the vector

      else if(translated[i] == '*'){

        //only adds to a valid string
        if(v.size() > 0 && v.back()[0] == 'M'){
          v.back() += translated[i];
          v.push_back("");
        }

        else{
          v.push_back("");
        }
      }

      // only adds the chars to a valid string
      else if(v.size() > 0 && v.back()[0] == 'M' && translated[i] != '*'){
        v.back() += translated[i];
      }
    }
  }

  int max = 0;
  int maxindex = 0;


  for(int i = 0; i < static_cast<int>(v.size()); i++){

    // goes over the vector and looks for the longest string in the vector
    if(static_cast<int>( v[i].size() ) > max && v[i].front() == 'M' && v[i].back() == '*'){

      //assigns max size of string
      // and index to a variable
      max = static_cast<int>(v[i].size());

      maxindex = i;

    }
  }

  //returns the largest string using the index of the largest string
  return v[maxindex];
}





#endif

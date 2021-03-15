// Copyright 2017 Martin C. Frith

#ifndef LAST_HH
#define LAST_HH

#include "Alphabet.hh"
#include "CyclicSubsetSeed.hh"
#include "MultiSequence.hh"
#include "SequenceFormat.hh"
#include "dna_words_finder.hh"
#include "qualityScoreUtil.hh"

namespace cbrc {

typedef MultiSequence::indexT indexT;

inline void err(const char *s) { throw std::runtime_error(s); }

// Read the next sequence, adding it to the MultiSequence
inline std::istream &appendSequence(MultiSequence &m, std::istream &in,
				    indexT maxSeqLen, sequenceFormat::Enum f,
				    const Alphabet &a, bool isKeepLowercase,
				    bool isMaskLowercase) {
  if (m.finishedSequences() == 0) maxSeqLen = -1;

  size_t oldSize = m.seqBeg(m.finishedSequences());

  if (f == sequenceFormat::fasta) {
    m.appendFromFasta(in, maxSeqLen);
  } else if (f == sequenceFormat::fastx) {
    m.appendFromFastx(in, maxSeqLen, false);
  } else if (f == sequenceFormat::fastxKeep) {
    m.appendFromFastx(in, maxSeqLen, true);
  } else if (f == sequenceFormat::prb) {
    m.appendFromPrb(in, maxSeqLen, a.size, a.decode);
  } else if (f == sequenceFormat::pssm) {
    m.appendFromPssm(in, maxSeqLen, a.encode, isMaskLowercase);
  } else {
    m.appendFromFastq(in, maxSeqLen, true);
  }

  if (!m.isFinished() && m.finishedSequences() == 0) {
    err("encountered a sequence that's too long");
  }

  size_t newSize = m.seqBeg(m.finishedSequences());

  // encode the newly-read sequence
  a.tr(m.seqWriter() + oldSize, m.seqWriter() + newSize, isKeepLowercase);

  if (isPhred(f)) {
    checkQualityCodes(m.qualityReader() + oldSize,
		      m.qualityReader() + newSize, qualityOffset(f));
  }  // assumes one quality code per letter

  return in;
}

inline void makeWordsFinder(DnaWordsFinder &wordsFinder,
			    const CyclicSubsetSeed *seeds, size_t numOfSeeds,
			    const uchar *lettersToNumbers,
			    bool isMaskLowercase) {
  wordsFinder.wordLength = 0;
  size_t wordLength = maxRestrictedSpan(seeds, numOfSeeds);
  if (wordLength) {
    if (numOfSeeds > dnaWordsFinderNull) {
      err("I can't handle so many word-restricted seed patterns, sorry");
    }
    if (wordLength >= CHAR_BIT * sizeof(unsigned) / 2) {
      err("I can't handle such long word restrictions in seed patterns, sorry");
    }
    assert(numOfSeeds > 0);
    size_t lengthOfAllWords = numOfSeeds * wordLength;
    std::vector<uchar> dnaMatches(lengthOfAllWords);
    for (size_t i = 0; i < numOfSeeds; ++i) {
      seeds[i].matchingDna(&dnaMatches[i * wordLength], wordLength);
    }
    wordsFinder.set(wordLength, numOfSeeds, &dnaMatches[0],
		    lettersToNumbers, isMaskLowercase);
  }
}

}

#endif

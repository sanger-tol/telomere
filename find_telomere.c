#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <zlib.h>

std::string rc(std::string str) {
  std::string DNAseq(str);

  size_t c = 0;
  // reverse
  for(int i = str.length()-1; i>=0; i--) {
     DNAseq[c++] = str[i];
  }
  // complement
  for (std::size_t i = 0; i < DNAseq.length(); ++i){
        switch (DNAseq[i]){
        case 'A':
            DNAseq[i] = 'T';
            break;
        case 'C':
            DNAseq[i] = 'G';
            break;
        case 'G':
            DNAseq[i] = 'C';
            break;
        case 'T':
            DNAseq[i] = 'A';
            break;
        }
    }
    return DNAseq;
}

static void find_motif(std::string str, std::string name, std::string str2) {
  std::size_t pos = str.find(str2);
  while (pos != std::string::npos) {
     std::size_t len = 0;
     std::cout << name << "\t" << str.length() << "\t0\t" << pos;
     while (str.substr(pos, str2.length()) == str2) {
        len+=str2.length();
        pos+=str2.length();
     }
     std::cout << "\t" << pos << "\t" << len<< std::endl;
     pos = str.find(str2, pos+1);
  }

  // now ref search
  std::string rev = rc(str2);
  pos = 0;
  pos = str.find(rev);
  while (pos != std::string::npos) {
     std::size_t len = 0;
     std::cout << name << "\t" << str.length() << "\t1\t" << pos;
     while (str.substr(pos, rev.length()) == rev) {
        len+=rev.length();
        pos+=rev.length();
     }
     std::cout << "\t" << pos << "\t" << len << std::endl;
     pos = str.find(rev, pos+1);
  }
}

static bool file_is_gzip(const char *path) {
  FILE *f = std::fopen(path, "rb");
  if (!f)
    return false;
  unsigned char b[2];
  std::size_t n = std::fread(b, 1, 2, f);
  std::fclose(f);
  return n == 2 && b[0] == 0x1f && b[1] == 0x8b;
}

static void trim_eol(std::string &line) {
  while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
    line.pop_back();
}

static void handle_line(const std::string &line, std::string &str, std::string &header,
                        int argc, char **argv) {
  if (line.find(">") == std::string::npos) {
    str.append(line);
  } else {
    if (str.length() > 0) {
      find_motif(str, header, (argc >= 3 ? argv[2] : "TTAGGG"));
    }
    str = "";
    header = line;
  }
}

int main(int argc, char * argv[])
{
  if (argc < 2) {
     std::cerr << "Error: invalid number of parameters" << std::endl;
     std::cerr << "Usage: find <input fasta> [optional sequence to search for, default is vertebrate TTAGG" << std::endl;
     exit(1);
   }

  std::string str ("");
  std::string header;

  if (file_is_gzip(argv[1])) {
    gzFile gz = gzopen(argv[1], "rb");
    if (!gz) {
      std::perror(argv[1]);
      return 1;
    }
    char buf[65536];
    while (gzgets(gz, buf, sizeof(buf)) != nullptr) {
      std::string line(buf);
      trim_eol(line);
      handle_line(line, str, header, argc, argv);
    }
    if (!gzeof(gz)) {
      int gzerrnum;
      const char *gzmsg = gzerror(gz, &gzerrnum);
      std::cerr << argv[1] << ": " << (gzmsg ? gzmsg : "gzip read error") << std::endl;
      gzclose(gz);
      return 1;
    }
    gzclose(gz);
  } else {
    std::ifstream infile(argv[1]);
    if (!infile) {
      std::perror(argv[1]);
      return 1;
    }
    std::string line;
    while (std::getline(infile, line)) {
      handle_line(line, str, header, argc, argv);
    }
    infile.close();
  }

  if (str.length() > 0) {
    find_motif(str, header, (argc >= 3 ? argv[2] : "TTAGGG"));
  }

  return 0;
}

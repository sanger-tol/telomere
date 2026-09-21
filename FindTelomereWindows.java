import java.io.BufferedReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;

public class FindTelomereWindows {
  private static final int WINDOW_SIZE = 1000;
  private static final int MIN_OFFSET = 0;
  private static final double DEFAULT_THRESHOLD = 0.4;
  private static final double DEFAULT_THRESHOLD_SPLIT = 0.05;
  private static double THRESHOLD = DEFAULT_THRESHOLD;

  private static final Comparator<Interval> INTERVAL_ORDER = new Comparator<Interval>() {
    public int compare(Interval a, Interval b) {
      int c = a.chrom.compareTo(b.chrom);
      if (c != 0) {
        return c;
      }
      if (a.start != b.start) {
        return Integer.compare(a.start, b.start);
      }
      if (a.end != b.end) {
        return Integer.compare(a.end, b.end);
      }
      return Integer.compare(a.strand, b.strand);
    }
  };

  private static class Interval {
    final String chrom;
    final int start;
    final int end;
    final int score;
    final int strand; // 0 = +, 1 = -, -1 = n/a (windows)

    Interval(String chrom, int start, int end, int score, int strand) {
      this.chrom = chrom;
      this.start = start;
      this.end = end;
      this.score = score;
      this.strand = strand;
    }
  }

  public FindTelomereWindows() {
  }

  public static void printUsage() {
    System.err.println("Usage: java FindTelomereWindows [--split] <in.telomere> <identity> [threshold]");
    System.err.println("No --split: writes <prefix>.telomere.bed and <prefix>.windows");
    System.err.println("--split: also writes <prefix>.fwd/.rev.telomere.bed and <prefix>.fwd/.rev.windows");
    System.err.println("All BED and windows outputs are sorted by chrom, start, end.");
    System.err.println("Windows are bedGraph-like: chrom start end score (score=ceil(1000*density), max 1000).");
    System.err.println("Telomere BED: chrom start end score strand.");
    System.err.println("");
  }

  public static void collectWindows(String name, BitSet b, int length, ArrayList<Interval> out) {
    if (b == null || out == null || name == null) {
      return;
    }

    for (int i = MIN_OFFSET; i <= length; i += WINDOW_SIZE / 5) {
      int car = b.get(i, Math.min(length, i + WINDOW_SIZE)).cardinality();
      int den = Math.min(WINDOW_SIZE, length - i);
      double density = (double) car / den;
      if (density >= THRESHOLD) {
        int score = (int) Math.ceil(density * 1000.0);
        if (score > 1000) {
          score = 1000;
        }
        out.add(new Interval(name, i, i + den, score, -1));
      }

      if (i + WINDOW_SIZE >= length) {
        break;
      }
    }
  }

  public static String outputPrefix(String inFile) {
    if (inFile.endsWith(".telomere")) {
      return inFile.substring(0, inFile.length() - ".telomere".length());
    }
    return inFile;
  }

  private static int bedScore(int start, int end) {
    int len = end - start;
    return len > 1000 ? 1000 : Math.max(0, len);
  }

  /**
   * Parse one find_telomere line into [chrom, length, strand, start, end].
   * Prefers tab-separated columns (name may contain spaces). Falls back to
   * taking the last 5 whitespace tokens as numerics and the first token as chrom.
   */
  private static String[] parseTelomereFields(String line) {
    line = line.trim();
    if (line.isEmpty()) {
      return null;
    }

    String chrom;
    String length;
    String strand;
    String start;
    String end;

    if (line.indexOf('\t') >= 0) {
      String[] tab = line.split("\t", -1);
      if (tab.length < 6) {
        return null;
      }
      int base = tab.length - 6;
      chrom = tab[base].trim().split("\\s+")[0];
      length = tab[base + 1].trim();
      strand = tab[base + 2].trim();
      start = tab[base + 3].trim();
      end = tab[base + 4].trim();
    } else {
      String[] ws = line.split("\\s+");
      if (ws.length < 6) {
        return null;
      }
      chrom = ws[0];
      length = ws[ws.length - 5];
      strand = ws[ws.length - 4];
      start = ws[ws.length - 3];
      end = ws[ws.length - 2];
    }

    return new String[] { chrom, length, strand, start, end };
  }

  private static void writeSortedBeds(ArrayList<Interval> beds, PrintStream out) {
    Collections.sort(beds, INTERVAL_ORDER);
    for (Interval iv : beds) {
      char strandChar = (iv.strand == 0) ? '+' : '-';
      out.format(Locale.US, "%s\t%d\t%d\t%d\t%c%n",
          iv.chrom, iv.start, iv.end, iv.score, strandChar);
    }
  }

  private static void writeSortedWindows(ArrayList<Interval> wins, PrintStream out) {
    Collections.sort(wins, INTERVAL_ORDER);
    for (Interval iv : wins) {
      out.format(Locale.US, "%s\t%d\t%d\t%d%n", iv.chrom, iv.start, iv.end, iv.score);
    }
  }

  public static void main(String[] args) throws Exception {
    boolean splitWindows = false;
    ArrayList<String> positionalArgs = new ArrayList<String>();
    for (int i = 0; i < args.length; i++) {
      String arg = args[i];
      if ("--split".equals(arg)) {
        splitWindows = true;
      } else {
        positionalArgs.add(arg);
      }
    }

    if (positionalArgs.size() < 2) {
      printUsage();
      System.exit(1);
    }

    THRESHOLD = splitWindows ? DEFAULT_THRESHOLD_SPLIT : DEFAULT_THRESHOLD;
    if (positionalArgs.size() >= 3) {
      THRESHOLD = Double.parseDouble(positionalArgs.get(2));
    }

    Double identity = Double.parseDouble(positionalArgs.get(1)) / 100;
    THRESHOLD = THRESHOLD * Math.pow(identity, 6);
    System.err.println("Given error rate of " + identity + " running with adjusted threshold of "
        + THRESHOLD + " due to survival prob " + Math.pow(identity, 6));

    String prefix = outputPrefix(positionalArgs.get(0));

    ArrayList<Interval> allBeds = new ArrayList<Interval>();
    ArrayList<Interval> fwdBeds = splitWindows ? new ArrayList<Interval>() : null;
    ArrayList<Interval> revBeds = splitWindows ? new ArrayList<Interval>() : null;
    ArrayList<Interval> allWins = new ArrayList<Interval>();
    ArrayList<Interval> fwdWins = splitWindows ? new ArrayList<Interval>() : null;
    ArrayList<Interval> revWins = splitWindows ? new ArrayList<Interval>() : null;

    if (splitWindows) {
      System.err.println("Writing sorted: " + prefix + ".telomere.bed, "
          + prefix + ".fwd.telomere.bed, " + prefix + ".rev.telomere.bed, "
          + prefix + ".windows, " + prefix + ".fwd.windows, " + prefix + ".rev.windows");
    } else {
      System.err.println("Writing sorted: " + prefix + ".telomere.bed, " + prefix + ".windows");
    }

    BufferedReader bf = Utils.getFile(positionalArgs.get(0), "telomere");
    BitSet fwd = null;
    BitSet rev = null;
    BitSet combined = null;
    String name = null;
    int length = 0;
    String line = null;

    while ((line = bf.readLine()) != null) {
      String[] fields = parseTelomereFields(line);
      if (fields == null) {
        continue;
      }

      String chrom;
      int seqLen;
      int strand;
      int start;
      int end;
      try {
        chrom = fields[0];
        seqLen = Integer.parseInt(fields[1]);
        strand = Integer.parseInt(fields[2]);
        start = Integer.parseInt(fields[3]);
        end = Integer.parseInt(fields[4]);
      } catch (NumberFormatException e) {
        System.err.println("Skipping malformed telomere line: " + line);
        continue;
      }

      if (combined == null || !chrom.equalsIgnoreCase(name)) {
        if (splitWindows) {
          collectWindows(name, fwd, length, fwdWins);
          collectWindows(name, rev, length, revWins);
        }
        collectWindows(name, combined, length, allWins);

        length = seqLen;
        name = chrom;
        combined = new BitSet(length);
        if (splitWindows) {
          fwd = new BitSet(length);
          rev = new BitSet(length);
        } else {
          fwd = null;
          rev = null;
        }
      }

      Interval bed = new Interval(chrom, start, end, bedScore(start, end), strand);
      allBeds.add(bed);
      if (splitWindows) {
        if (strand == 0) {
          fwdBeds.add(bed);
          fwd.set(start, end);
        } else if (strand == 1) {
          revBeds.add(bed);
          rev.set(start, end);
        } else {
          System.err.println("Unknown orientation on line " + line);
          System.exit(1);
        }
      } else if (strand != 0 && strand != 1) {
        System.err.println("Unknown orientation on line " + line);
        System.exit(1);
      }
      combined.set(start, end);
    }

    if (splitWindows) {
      collectWindows(name, fwd, length, fwdWins);
      collectWindows(name, rev, length, revWins);
    }
    collectWindows(name, combined, length, allWins);
    bf.close();

    PrintStream allBedOut = new PrintStream(prefix + ".telomere.bed");
    writeSortedBeds(allBeds, allBedOut);
    allBedOut.close();

    PrintStream allWinOut = new PrintStream(prefix + ".windows");
    writeSortedWindows(allWins, allWinOut);
    allWinOut.close();

    if (splitWindows) {
      PrintStream fwdBedOut = new PrintStream(prefix + ".fwd.telomere.bed");
      writeSortedBeds(fwdBeds, fwdBedOut);
      fwdBedOut.close();

      PrintStream revBedOut = new PrintStream(prefix + ".rev.telomere.bed");
      writeSortedBeds(revBeds, revBedOut);
      revBedOut.close();

      PrintStream fwdWinOut = new PrintStream(prefix + ".fwd.windows");
      writeSortedWindows(fwdWins, fwdWinOut);
      fwdWinOut.close();

      PrintStream revWinOut = new PrintStream(prefix + ".rev.windows");
      writeSortedWindows(revWins, revWinOut);
      revWinOut.close();
    }
  }
}

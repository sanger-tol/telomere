import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.io.File;
import java.io.PrintStream;
import java.util.BitSet;

public class FindTelomereWindows {  
  private static final NumberFormat nf = new DecimalFormat("############.#");
  private static final int WINDOW_SIZE = 1000;
  private static final int MIN_OFFSET = 0;
  private static final double DEFAULT_THRESHOLD = 0.4;
  private static final double DEFAULT_THRESHOLD_SPLIT = 0.05;
  private static double THRESHOLD = DEFAULT_THRESHOLD;

   public FindTelomereWindows() {
   }

   public static void printUsage() {
      System.err.println("Usage: java -jar FindTelomereWindows.jar [--split] <in> <identity> [threshold]");
      System.err.println("This program sizes a fasta or fastq file. Multiple fasta files can be supplied by using a comma-separated list.");
      System.err.println("Example usage: getHist fasta1.fasta,fasta2.fasta");
      System.err.println("");
   }
  public static void processScaffoldToStream(String name, BitSet b, int length, PrintStream out) {
      if (b == null) { return; }

      for (int i = MIN_OFFSET; i <= length; i+=WINDOW_SIZE/5) {
         int car = b.get(i, Math.min(length, i+WINDOW_SIZE)).cardinality();
         int den = Math.min(WINDOW_SIZE, length-i);
         if ((double)car / den >= THRESHOLD)
            out.println("Window\t" + name + "\t" + length + "\t" + i + "\t" + (i+den) + "\t" + ((double)car / den));

         if (i+WINDOW_SIZE >= length)
            break;
      }
  }
   public static void processScaffold(String name, BitSet b, int length) {
      processScaffoldToStream(name, b, length, System.out);
   }
   
  public static String splitOutputPrefix(String inFile) {
      if (inFile.endsWith(".telomere")) {
         return inFile.substring(0, inFile.length() - ".telomere".length());
      }
      return inFile;
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

      if (positionalArgs.size() < 2) { printUsage(); System.exit(1);}

      THRESHOLD = splitWindows ? DEFAULT_THRESHOLD_SPLIT : DEFAULT_THRESHOLD;
      if (positionalArgs.size() >= 3) {
          THRESHOLD = Double.parseDouble(positionalArgs.get(2));
      }

      Double identity = Double.parseDouble(positionalArgs.get(1)) / 100;
      THRESHOLD = THRESHOLD * Math.pow(identity, 6);
      System.err.println("Given error rate of " + identity + " running with adjusted threshold of " + THRESHOLD + " due to survival prob " + Math.pow(identity, 6));

      BufferedReader bf = Utils.getFile(positionalArgs.get(0), "telomere");
      if (splitWindows) {
         String splitPrefix = splitOutputPrefix(positionalArgs.get(0));
         PrintStream fwdOut = new PrintStream(splitPrefix + ".fwd.windows");
         PrintStream revOut = new PrintStream(splitPrefix + ".rev.windows");
         System.err.println("Writing split windows to " + splitPrefix + ".fwd.windows and " + splitPrefix + ".rev.windows");
         BitSet fwd = null;
         BitSet rev = null;
         String name = null;
         int length = 0;
         int offset = 0;
         String line = null;
         while ((line = bf.readLine()) != null) {
             String[] split = line.trim().split("\\s+");
             if (split.length > 6) { offset = split.length - 6; }
             if (fwd == null || !split[0].equalsIgnoreCase(name)) {
                processScaffoldToStream(name, fwd, length, fwdOut);
                processScaffoldToStream(name, rev, length, revOut);
                fwd = new BitSet(Integer.parseInt(split[1+offset]));
                rev = new BitSet(Integer.parseInt(split[1+offset]));
                length = Integer.parseInt(split[1+offset]);
                name = split[0];
             }
             if (Integer.parseInt(split[2+offset]) == 0) {
                fwd.set(Integer.parseInt(split[3+offset]), Integer.parseInt(split[4+offset]));
             } else if (Integer.parseInt(split[2+offset]) == 1) {
                rev.set(Integer.parseInt(split[3+offset]), Integer.parseInt(split[4+offset]));
             } else {
                System.err.println("Unknown orientation on line " + line);
                System.exit(1);
             }
          }
         processScaffoldToStream(name, fwd, length, fwdOut);
         processScaffoldToStream(name, rev, length, revOut);
         fwdOut.close();
         revOut.close();
      } else {
         // initialize sizes
         BitSet scaffold = null;
         String name = null;
         int length = 0;
         String line = null;
         while ((line = bf.readLine()) != null) {
             String[] split = line.trim().split("\\s+");
             if (scaffold == null || !split[0].equalsIgnoreCase(name)) {
                processScaffold(name, scaffold, length);
                length = Integer.parseInt(split[split.length-5]);
                scaffold = new BitSet(length);
                name = split[0];
             }
             // ignoring strandedness in default mode
             scaffold.set(Integer.parseInt(split[split.length-3]), Integer.parseInt(split[split.length-2]));
         }
         processScaffold(name, scaffold, length);
      }
      bf.close();
   }
}

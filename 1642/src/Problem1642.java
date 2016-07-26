import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StreamTokenizer;

public class Problem1642 {

	private static final int INF = 100000;

	public static void main(String[] args) {
		try {
			run();
		} catch (Exception e) {
			System.err.println(e);
		}
	}
	
	private static void run() throws FileNotFoundException, IOException {
		boolean oj = System.getProperty("ONLINE_JUDGE") != null;
		Reader reader = oj ? new InputStreamReader(System.in) : new FileReader("input.txt");
		// Writer writer = oj ? new OutputStreamWriter(System.out) : new FileWriter("output.txt");
		StreamTokenizer in = new StreamTokenizer(new BufferedReader(reader));
		// PrintWriter out = new PrintWriter(writer);
		
		int n = nextInt(in);
		int exit = nextInt(in);
		int[] numbers = new int[n];
		for (int i = 0; i < n; ++i) {
			numbers[i] = nextInt(in);
		}
		
		int numForward = getNumMoves(exit, numbers, 1);
		if (numForward != INF) {
			int numBackward = getNumMoves(exit, numbers, -1);
			System.out.print(numForward);
			System.out.print(" ");
			System.out.print(numBackward);
		} else {
			System.out.print("Impossible");
		}
	}

	private static int getNumMoves(int exit, int[] numbers, int inc) {
		int now = 0;
		int moves = 0;
		while (now != exit && moves < INF) {
			++moves;
			now += inc;
			int j = 0;
			while (j < numbers.length && now != numbers[j]) {
				++j;
			}
			if (j < numbers.length) {
				inc = -inc;
			}
		}
		return moves;
	}

	private static int nextInt(StreamTokenizer in) throws IOException {
		in.nextToken();
		return (int) in.nval;
	}
}

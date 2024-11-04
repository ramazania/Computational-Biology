/**
  This program will take in a variable number of command line arguments
  and will print the total number of arguments, followed by each one
  to the screen.  Here is an example of how to run it and what the
  output might look like:

  $ javac Echo.java
  $ java Echo 1stArg 2ndArg 3rdArg
    Total number of arguments: 3
    Argument 0 : 1stArg
    Argument 1 : 2ndArg
    Argument 2 : 3rdArg
 **/

public class Echo {

    public static void main(String[] args) 
    {
            int numArgs = args.length;
            
            System.out.println("Total number of arguments: " + numArgs);

            for (int i = 0; i < numArgs; i++)
            {
                    System.out.println("Argument " + i + " : " + args[i]);
            }
    
    }

}

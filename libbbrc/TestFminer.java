class TestFminer {
    public static void main(String[] args) {
       System.loadLibrary("fminer");
       Fminer MyFminer;
       MyFminer = new Fminer();
       MyFminer.AddCompound ("COC1=CC=C(C=C1)C2=NC(=C([NH]2)C3=CC=CC=C3)C4=CC=CC=C4", 1);
       MyFminer.AddCompound ("O=C1NC(=S)NC(=O)C1C(=O)NC2=CC=CC=C2", 2);
          // ... continue adding compounds
       MyFminer.AddActivity((boolean) true, 1);
       MyFminer.AddActivity((boolean) false, 2);
          // ... continue adding activities (true for active, false for inactive)
       System.out.println(MyFminer.GetNoCompounds() + " compounds");
       // Toy example: special settings for mining all fragments
       MyFminer.SetChisqSig(0); // use no significance constraint
       MyFminer.SetRefineSingles(true); // refine structures with support 1
       // gather results for every root node in vector instead of immediate output
       MyFminer.SetConsoleOut(false);
       for (int j = 0; j < (int) MyFminer.GetNoRootNodes(); j++)
       {
          SVector result = MyFminer.MineRoot(j);
          for(int i = 0; i < result.size(); i++)
          {
            System.out.println(result.get(i));
          }
       }
       MyFminer = null;
    }
}

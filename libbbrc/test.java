public class test {
    public static void main(String args[]) {
       System.loadLibrary("bbrc");
       Bbrc MyFminer;
       MyFminer = new Bbrc();
       // Toy example: special settings for mining all fragments
       MyFminer.SetChisqSig(0); // use no significance constraint
       MyFminer.SetRefineSingles(true); // refine structures with support 1
       MyFminer.SetConsoleOut(false);
       // Add compounds below. IMPORTANT! DO NOT CHANGE SETTINGS AFTER ADDING COMPOUNDS!
       MyFminer.AddCompound ("COC1=CC=C(C=C1)C2=NC(=C([NH]2)C3=CC=CC=C3)C4=CC=CC=C4", 1);
       MyFminer.AddCompound ("O=C1NC(=S)NC(=O)C1C(=O)NC2=CC=CC=C2", 2);
       // ... continue adding compounds
       MyFminer.AddActivity(1.0F, 1);
       MyFminer.AddActivity(0.0F, 2);
       // ... continue adding activities (1.0F for active, 0.0F for inactive)
       System.out.println(MyFminer.GetNoCompounds() + " compounds");
       // gather results for every root node in vector instead of immediate output
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
 



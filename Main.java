import java.util.Scanner;      
  
public class Main {      
    public static void main(String[] args) {      
        Scanner scanner = new Scanner(System.in);      
        int i = scanner.nextInt();               
        int reversed = 0;      
            
        while (i != 0) {      
            int digit = i % 10;     
            reversed = reversed * 10 + digit;     
            i /= 10;    
        }    
        System.out.println(reversed);      
        scanner.close(); 
    }      
}
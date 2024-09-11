import java.util.Scanner;  
  
public class Time{  
    public static void main(String[] args) {  
        Scanner scanner = new Scanner(System.in);  
    
        int firstHour = scanner.nextInt();  
        int firstMinute = scanner.nextInt();  
   
        int secondHour = scanner.nextInt();  
        int secondMinute = scanner.nextInt();  
  
        int diffHours = secondHour - firstHour;  
        int diffMinutes = secondMinute - firstMinute;  
  
        if (diffMinutes < 0) {  
            diffMinutes += 60;  
            diffHours--;  
        }  
  
        System.out.println(diffHours+" "+diffMinutes);  
  
        scanner.close();  
    }  
}
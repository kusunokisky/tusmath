package sjtools;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.apache.poi.xssf.usermodel.XSSFCell;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

class DataConvertExcel{
	private XSSFWorkbook wb;
	private XSSFSheet sheet;
	private XSSFRow row;
	private XSSFCell cell;
	private FileOutputStream out = null;

	DataConvertExcel() {
		wb = new XSSFWorkbook();
		sheet = wb.createSheet();
	}
	void setData(double data,int nrow,int ncolumn){
		row = sheet.getRow(nrow);
		if(row == null){
			row = sheet.createRow(nrow);
			cell = row.createCell(ncolumn);
		}else{
			cell = row.getCell(ncolumn);
			if(cell == null){
				cell = row.createCell(ncolumn);
			}
		}
		cell.setCellValue(data);
	}
	/**
	 * データを任意の列に対して書き込んでいきます
	 * @param data 書き込むデータ
	 * @param nrow 書き込む行
	 * @param ncolumn 書き込む列
	 */
	void setData(double[] data,int nrow,int ncolumn){
		for(int i = 0;i < data.length;i++,nrow++){
			setData(data[i], nrow, ncolumn);
		}
	}
	/**
	 * 現在のデータを基に.xlsx形式でファイルを保存します
	 */
	void createDataFile(){
		JFileChooser filechooser = new JFileChooser();
		int selected = filechooser.showSaveDialog(null);
		if(selected == JFileChooser.APPROVE_OPTION){
			File file = filechooser.getSelectedFile();
			try{
				out = new FileOutputStream(file.getAbsolutePath());
				wb.write(out);
			}catch(IOException ioe){
				System.out.println(ioe.toString());
			}
			finally{
				try{
					out.close();
				}catch(IOException ioe_){
					System.out.println(ioe_.toString());
				}
			}
		}else if(selected == JFileChooser.CANCEL_OPTION){
			JOptionPane.showMessageDialog(null, "作業がキャンセルされました");
		}


	}
}

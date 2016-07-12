package gui;

import java.awt.BorderLayout;
import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.JLabel;
import java.awt.Font;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.JFileChooser;

import java.awt.event.ActionListener;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class Janela extends JFrame {

	private JPanel contentPane;
	private JTextField textField;
	private JTextField txtCamera;
	private JTextField txtObjeto;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					Janela frame = new Janela();
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
	public Janela() {
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		setBounds(100, 100, 358, 297);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		this.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				try{
					PrintWriter arquivo = new PrintWriter("config.txt", "UTF-8");
					arquivo.println("q");
					arquivo.close();
					Janela.this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				} catch (FileNotFoundException | UnsupportedEncodingException e1) {
					e1.printStackTrace();
				}
			}
		});
		
		JLabel lblnguloDeRotao = new JLabel("\u00C2ngulo de Rota\u00E7\u00E3o");
		lblnguloDeRotao.setBounds(10, 24, 115, 14);
		contentPane.add(lblnguloDeRotao);
		
		textField = new JTextField();
		textField.setBounds(169, 21, 115, 20);
		contentPane.add(textField);
		textField.setColumns(10);
		
		JLabel lblRotao = new JLabel("Rota\u00E7\u00E3o");
		lblRotao.setBounds(149, 64, 71, 14);
		contentPane.add(lblRotao);
		
		JButton btnOz = new JButton("Eixo OZ");
		btnOz.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				fazer_rotacao("oz");
			}
		});
		btnOz.setBounds(10, 97, 89, 23);
		contentPane.add(btnOz);
		
		JButton btnEixoOy = new JButton("Eixo OY");
		btnEixoOy.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				fazer_rotacao("oy");
			}
		});
		btnEixoOy.setBounds(131, 97, 89, 23);
		contentPane.add(btnEixoOy);
		
		JButton btnEixoOz = new JButton("Eixo OX");
		btnEixoOz.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				fazer_rotacao("ox");
			}
		});
		btnEixoOz.setBounds(247, 97, 89, 23);
		contentPane.add(btnEixoOz);
		
		JLabel lblMudanaDeObjeto = new JLabel("Mudan\u00E7a de Objeto");
		lblMudanaDeObjeto.setBounds(120, 131, 115, 14);
		contentPane.add(lblMudanaDeObjeto);
		
		JLabel lblCmera = new JLabel("C\u00E2mera");
		lblCmera.setBounds(10, 163, 46, 14);
		contentPane.add(lblCmera);
		
		txtCamera = new JTextField();
		txtCamera.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent arg0) {
				JFileChooser fc = new JFileChooser();
				FileNameExtensionFilter filtro = new FileNameExtensionFilter("Arquivos CFG", "cfg");
				fc.setFileFilter(filtro);
				int val = fc.showOpenDialog(Janela.this);
				if(val == JFileChooser.APPROVE_OPTION){
					String caminho = fc.getSelectedFile().getAbsolutePath();
					txtCamera.setText(caminho);
				}
			}
		});
		txtCamera.setBounds(66, 160, 270, 20);
		contentPane.add(txtCamera);
		txtCamera.setColumns(10);
		
		JLabel lblObjeto = new JLabel("Objeto");
		lblObjeto.setBounds(10, 190, 46, 14);
		contentPane.add(lblObjeto);
		
		txtObjeto = new JTextField();
		txtObjeto.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent arg0) {
				JFileChooser fc = new JFileChooser();
				FileNameExtensionFilter filtro = new FileNameExtensionFilter("Arquivos BYU", "byu");
				fc.setFileFilter(filtro);
				int val = fc.showOpenDialog(Janela.this);
				if(val == JFileChooser.APPROVE_OPTION){
					String caminho = fc.getSelectedFile().getAbsolutePath();
					txtObjeto.setText(caminho);
				}
			}
		});
		txtObjeto.setColumns(10);
		txtObjeto.setBounds(66, 187, 270, 20);
		contentPane.add(txtObjeto);
		
		JButton btnMudarObjeto = new JButton("Mudar Objeto");
		btnMudarObjeto.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				mudar_objeto();
			}
		});
		btnMudarObjeto.setBounds(120, 225, 115, 23);
		contentPane.add(btnMudarObjeto);
	}
	
	public void fazer_rotacao(String eixo){
		String str_angulo = this.textField.getText();
		//float angulo = Float.parseFloat(str_angulo);
		try {
			PrintWriter arquivo = new PrintWriter("config.txt", "UTF-8");
			arquivo.println("r");
			arquivo.println(str_angulo);
			arquivo.println(eixo);
			arquivo.close();
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		this.dispose();
	}
	
	public void mudar_objeto(){
		String cam_camera = this.txtCamera.getText();
		String cam_objeto = this.txtObjeto.getText();
		try {
			PrintWriter arquivo = new PrintWriter("config.txt", "UTF-8");
			arquivo.println("m");
			arquivo.println(cam_camera);
			arquivo.println(cam_objeto);
			arquivo.close();
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		this.dispose();
	}
}

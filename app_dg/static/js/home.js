const buttonClick = () => {
	const job_form = document.getElementById("theForm")
	const job_button = document.getElementById("job_button")

	job_button.classList.remove('opacity-100');
	job_button.classList.add('opacity-0');
	job_button.classList.add('hidden');

	job_form.classList.remove("invisible");
	job_form.classList.remove("absolute", "bottom-0", "left-0");
	job_form.classList.add('opacity-100');
}

const theFile = document.getElementById("theFile");
const checkMail = (event) => {
	const res = String(event.target.value).toLowerCase().match(
	  /^(([^<>()[\]\\.,;:\s@"]+(\.[^<>()[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/
	);
	const button_label = document.getElementById("button_label");

	if (res == null) {
	  event.target.classList.remove("text-green-600");
	  event.target.classList.add("text-red-500");
	  event.target.valid = false;
	  button_label.classList.remove("hover:bg-green-600","hover:text-white", "text-green-600", "cursor-pointer")
	  button_label.classList.add("bg-gray-600","text-white")
	  theFile.disabled = true;

	} else {
	  event.target.classList.remove("text-red-500");
	  event.target.classList.add("text-green-600");
	  event.target.valid = true;
	  button_label.classList.remove("bg-gray-600","text-white")
	  button_label.classList.add("hover:bg-green-600","hover:text-white", "text-green-600", "cursor-pointer")
	  theFile.disabled = false;
	}

}

const email = document.getElementById("theMail");
email.value = ""
email.addEventListener('input', checkMail);


const sendRequest = (event) => {
	let request = new XMLHttpRequest();
	request.open('POST', '');
	request.send(formdata);
	request.onreadystatechange = () => {
		if (request.readyState == XMLHttpRequest.DONE) {
			var OK = 200;

			if (request.status === OK) {
				window.location.href = request.responseURL;
			}
			else {
				console.log ('Error: ' + request.status); 
			}
		}
	};
}

theFile.addEventListener('change', sendRequest);

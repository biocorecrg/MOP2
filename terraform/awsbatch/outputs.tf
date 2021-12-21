output "public_dns" {
  value = aws_instance.classroom.*.public_dns
}

output "instance_id" {
  value = aws_instance.classroom.*.id
}

output "bucket_name" {
  value = aws_s3_bucket.class-bucket.*.bucket
}

